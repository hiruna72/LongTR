// ShardPlanner -- predict, per catalog locus, how much compressed BAM data LongTR's region query
// will decompress, so the catalog can be sharded into pieces of equal predicted WORK rather than
// equal locus COUNT.
//
// WHY THIS EXISTS
// Shards cut by line count are badly unbalanced: cost per locus is driven by how far htslib has to
// backtrack to find the earliest read that could overlap, which depends on local read length and
// depth, not on how many loci you asked for. Two 94-locus shards measured on the chr21 fixture
// differed 17x in genomic span and 1.5x in runtime. With GNU parallel dispatching many shards, wall
// clock is set by the slowest one, so flattening that distribution buys real time. It does NOT
// reduce total CPU -- the same work is done, just spread evenly.
//
// WHY IT MUST LINK HTSLIB
// The prediction is the set of file ranges htslib will read, i.e. the chunk list that
// hts_itr_query() builds in iter->off[]. Reconstructing that outside htslib means reimplementing
// bin selection, linear-index pruning and chunk merging, and two attempts at that have now been
// measured against instrumented ground truth and failed badly (a hand-written BAI parser predicted
// 227 MB for a locus needing ~3,200 records; a linear-index-only proxy scored Spearman rho = 0.07
// against measured per-locus cost). Asking htslib itself is exact by construction.
//
// COST
// sam_itr_queryi() consults only the loaded index -- it reads no alignment records -- so planning is
// essentially free (~0.02 s for 200 loci) even though it predicts work that takes minutes to do.
// Only the BAM header and the .bai are read.
//
// WHAT THE NUMBER IS
// Compressed bytes spanned by the merged chunk list. Virtual offsets pack the BGZF block address in
// the high 48 bits, so (v>>16)-(u>>16) is the compressed distance a chunk covers. Measured against
// instrumented htslib this runs 2-6% high in absolute terms (it cannot see hts_itr_next()'s early
// termination) but tracks RELATIVE cost -- the only thing sharding needs -- to within ~3%.
//
// Intervals are merged per FILE: offsets from different BAMs live in different address spaces and
// must never be merged with each other, only summed.

#include <getopt.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "error.h"

static const int DEFAULT_PAD = 1000;   // must match BamProcessor::MAX_MATE_DIST

static void print_usage(){
  std::cerr << "Usage: ShardPlanner --bams <list> --regions <catalog.bed> [--pad <bp>]\n\n"
	    << "\t" << "--bams      <a.bam,b.bam>   " << "\t" << "Comma separated list of BAMs. Each needs its .bai alongside\n"
	    << "\t" << "--bam-files <bams.txt>      " << "\t" << "File listing BAMs, one per line. Alternative to --bams\n"
	    << "\t" << "--regions   <catalog.bed>   " << "\t" << "Catalog BED. Pass the SAME file LongTR will be given --\n"
	    << "\t" << "                            " << "\t" << " if the caller reorders the catalog, plan the reordered copy\n"
	    << "\t" << "--pad       <bp>            " << "\t" << "Query padding either side of each locus (Default = " << DEFAULT_PAD << ",\n"
	    << "\t" << "                            " << "\t" << " matching LongTR's --max-mate-dist)\n\n"
	    << "Writes one line per catalog locus to stdout, input order preserved:\n"
	    << "\t<predicted_compressed_bytes>\\t<original BED line>\n"
	    << "preceded by a '#' header line carrying the totals.\n\n"
	    << "The prediction is an upper bound in absolute terms; it is intended for comparing loci\n"
	    << "against each other, not for estimating absolute runtime.\n" << std::endl;
}

static std::vector<std::string> split_commas(const std::string& s){
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ','))
    if (!tok.empty()) out.push_back(tok);
  return out;
}

struct BamHandle {
  samFile*   fp;
  sam_hdr_t* hdr;
  hts_idx_t* idx;
};

int main(int argc, char** argv){
  std::string bams_string, bamfile_string, region_file;
  int pad = DEFAULT_PAD;

  static struct option long_options[] = {
    {"bams",      required_argument, 0, 'b'},
    {"bam-files", required_argument, 0, 'B'},
    {"regions",   required_argument, 0, 'r'},
    {"pad",       required_argument, 0, 'p'},
    {"help",      no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };
  while (true){
    int option_index = 0;
    int c = getopt_long(argc, argv, "b:B:r:p:h", long_options, &option_index);
    if (c == -1) break;
    switch(c){
    case 'b': bams_string    = std::string(optarg); break;
    case 'B': bamfile_string = std::string(optarg); break;
    case 'r': region_file    = std::string(optarg); break;
    case 'p': pad            = atoi(optarg);        break;
    case 'h': print_usage(); return 0;
    case '?': print_usage(); return 1;
    default:  break;
    }
  }
  if (argc == 1){ print_usage(); return 1; }

  std::vector<std::string> bam_paths = split_commas(bams_string);
  if (!bamfile_string.empty()){
    std::ifstream input(bamfile_string.c_str());
    if (!input.is_open())
      printErrorAndDie("Failed to open --bam-files file " + bamfile_string);
    std::string line;
    while (std::getline(input, line)){
      if (!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
      if (!line.empty()) bam_paths.push_back(line);
    }
  }
  if (bam_paths.empty()) printErrorAndDie("Must specify --bams or --bam-files");
  if (region_file.empty()) printErrorAndDie("Must specify --regions");
  if (pad < 0) printErrorAndDie("--pad cannot be negative");

  std::vector<BamHandle> bams;
  for (size_t i = 0; i < bam_paths.size(); i++){
    BamHandle h;
    h.fp = sam_open(bam_paths[i].c_str(), "r");
    if (h.fp == NULL) printErrorAndDie("Failed to open " + bam_paths[i]);
    h.hdr = sam_hdr_read(h.fp);
    if (h.hdr == NULL) printErrorAndDie("Failed to read the header for " + bam_paths[i]);
    // Index only -- no alignment data is read from the file at any point below.
    h.idx = sam_index_load(h.fp, bam_paths[i].c_str());
    if (h.idx == NULL) printErrorAndDie("Failed to load the index for " + bam_paths[i]);
    bams.push_back(h);
  }

  std::ifstream regions(region_file.c_str());
  if (!regions.is_open()) printErrorAndDie("Failed to open --regions file " + region_file);

  std::vector<std::string> out_lines;
  std::vector<uint64_t>    out_costs;
  uint64_t total = 0;
  int64_t  n_loci = 0, n_missing_contig = 0;
  std::string line;
  std::vector<std::pair<uint64_t,uint64_t> > iv;

  while (std::getline(regions, line)){
    if (!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
    if (line.empty() || line[0] == '#') continue;

    std::stringstream ss(line);
    std::string chrom; int64_t start, stop;
    if (!(ss >> chrom >> start >> stop))
      printErrorAndDie("Malformed BED line in " + region_file + ": " + line);

    // Exactly the window BamProcessor asks htslib for (bam_processor.cpp: region start/stop +/- MAX_MATE_DIST)
    int64_t beg = start - pad; if (beg < 0) beg = 0;
    int64_t end = stop + pad;

    uint64_t cost = 0;
    bool seen_contig = false;
    for (size_t i = 0; i < bams.size(); i++){
      int tid = sam_hdr_name2tid(bams[i].hdr, chrom.c_str());
      if (tid < 0) continue;                       // contig absent from this BAM
      seen_contig = true;
      hts_itr_t* itr = sam_itr_queryi(bams[i].idx, tid, beg, end);
      if (itr == NULL) continue;                   // no data for this region

      iv.clear();
      for (int k = 0; k < itr->n_off; k++)
	iv.push_back(std::make_pair((uint64_t)itr->off[k].u >> 16, (uint64_t)itr->off[k].v >> 16));
      hts_itr_destroy(itr);

      // Merge within this file only, then add. htslib walks chunks in offset order and never
      // revisits, so overlapping chunks would otherwise be counted twice.
      std::sort(iv.begin(), iv.end());
      for (size_t k = 0; k < iv.size(); ){
	uint64_t lo = iv[k].first, hi = iv[k].second;
	size_t j = k + 1;
	while (j < iv.size() && iv[j].first <= hi){
	  if (iv[j].second > hi) hi = iv[j].second;
	  j++;
	}
	if (hi > lo) cost += hi - lo;
	k = j;
      }
    }
    if (!seen_contig) n_missing_contig++;

    out_costs.push_back(cost);
    out_lines.push_back(line);
    total += cost;
    n_loci++;
  }

  if (n_missing_contig > 0)
    std::cerr << "WARNING: " << n_missing_contig << " catalog line(s) name a contig absent from every BAM; "
	      << "they are emitted with cost 0 and will be sharded but produce no calls" << std::endl;

  std::printf("#ShardPlanner\tpad=%d\tbams=%zu\tloci=%" PRId64 "\ttotal_bytes=%" PRIu64 "\n",
	      pad, bam_paths.size(), n_loci, total);
  for (size_t i = 0; i < out_lines.size(); i++)
    std::printf("%" PRIu64 "\t%s\n", out_costs[i], out_lines[i].c_str());

  for (size_t i = 0; i < bams.size(); i++){
    hts_idx_destroy(bams[i].idx);
    sam_hdr_destroy(bams[i].hdr);
    sam_close(bams[i].fp);
  }
  return 0;
}
