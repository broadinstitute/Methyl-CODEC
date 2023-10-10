//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "spoa/spoa.hpp"

#include "Alignment.h"
#include "AlignmentConsensus.h"
#include "BamIO.h"
#include <seqan/align.h>
#include "seqan_init_score.h"
#include "SeqLib/BWAWrapper.h"
#include "BamRecordExt.h"

extern "C" {
#include "bwa/ksw.h"
#ifdef USE_MALLOC_WRAPPERS
#  include "bwa/malloc_wrap.h"
#endif
}

#ifdef BGZF_MAX_BLOCK_SIZE
#pragma push_macro("BGZF_MAX_BLOCK_SIZE")
#undef BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE_BAK
#endif

#ifdef BGZF_BLOCK_SIZE
#pragma push_macro("BGZF_BLOCK_SIZE")
#undef BGZF_BLOCK_SIZE
#define BGZF_BLOCK_SIZE_BAK
#endif

#include "InsertSeqFactory.h"

#ifdef BGZF_MAX_BLOCK_SIZE_BAK
#undef BGZF_MAX_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_MAX_BLOCK_SIZE")
#endif

#ifdef BGZF_BLOCK_SIZE_BAK
#undef BGZF_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_BLOCK_SIZE")
#endif


using std::string;
struct CssOptions {
  string bam;
  int mapq = 10;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool allow_nonoverlapping_pair = false;
  bool clip3 = false;
//  int consensus_mode = 0;
  int pair_min_overlap = 1;
  int max_read = 0;
  int min_overlap = 30;
  bool trim_overhang = false;
  string reference = "";
  string tmpdir = "/tmp";
  //bool output_singleend = false;
  int minbq = 0;
  string outprefix = "";
  int thread = 1;
};


static struct option  consensus_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"reference",                required_argument,      0,        'r'},
    {"load_supplementary",       no_argument,            0,        'l'},
    {"clip3",                    no_argument,            0,        'C'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"baseq",                    required_argument ,     0,        'q'},
    {"outprefix",                required_argument ,     0,        'o'},
    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"trim_overhang",            no_argument,            0,        't'},
    {"allow_nonoverlapping_pair",no_argument,            0,        'i'},
    //{"output_singleend",         no_argument,            0,        's'},
    //hidden parameter
   // {"consensus_mode",           required_argument ,     0,        'M'},
    {"dirtmp",                   required_argument ,     0,        'd'},
    {"thread",                   required_argument,      0,        'T'},
    {"max_read",                   required_argument,      0,        'M'},
    {0,0,0,0}
};

const char* consensus_short_options = "b:m:M:o:lCp:q:d:tT:ir:";

void consensus_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              Input bam [required]\n";
  std::cerr<< "-o/--outprefix,                        Output sample prefix [required].\n";
  std::cerr<< "-m/--mapq,                             Min mapping quality [10].\n";
  std::cerr<< "-q/--baseq,                            If one of the baseq < cutoff, make all baseq = 2 so that VC will ingnore them. [0].\n";
  std::cerr<< "-r/--reference,                        Reference for alignment. [null].\n";
  //Supplementary not support currently
  //std::cerr<< "-l/--load_supplementary,               Include supplementary alignment [false].\n";
  std::cerr<< "-t/--trim_overhang,                    When perform paired-end consensus, if true then only do consensus of the overlapped region [false].\n";
  std::cerr<< "-C/--clip3,                            trim the 3'end soft clipping [false].\n";
//  std::cerr<< "-s/--output_singleend,                 The R1R2 consensus will be output in a single end format [false].\n";
  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair [1].\n";
  std::cerr<< "-d/--dirtmp,                           Temporary dir for sorted bam [/tmp]\n";
  std::cerr<< "-T/--thread,                           Number of threads for sort [1]\n";
  std::cerr<< "-i/--allow_nonoverlapping_pair,        Allow output of non-overlaping pairs, usually caused by intermolecular ligation. This will simply print the original reads.  [false]\n";
  std::cerr<< "-M/--max_read,                         Maximum number of read pair process. [Inf]\n";
}

int consensus_parse_options(int argc, char* argv[], CssOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, consensus_short_options, consensus_long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'b':
        opt.bam = optarg;
        break;
      case 'o':
        opt.outprefix = optarg;
        break;
      case 'r':
        opt.reference = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.minbq = atoi(optarg);
        break;
      case 'l':
        opt.load_supplementary = true;
        break;
      case 't':
        opt.trim_overhang = true;
        break;
      case 'i':
        opt.allow_nonoverlapping_pair = true;
        break;
      case 'C':
        opt.clip3 = true;
        break;
//      case 's':
//        opt.output_singleend = true;
//        break;
      case 'd':
        opt.tmpdir = optarg;
        break;
      case 'M':
        opt.max_read = atoi(optarg);
        break;
      case 'p':
        opt.pair_min_overlap = atoi(optarg);
        break;
      case 'T':
        opt.thread = atoi(optarg);
        break;
      default:consensus_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


int codec_ms_align(int argc, char ** argv) {
  CssOptions opt;
  int parse_ret =  consensus_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    consensus_print_help();
    return 1;
  }
  if (opt.outprefix.empty()) {
    std::cerr << "-o/--outprefix is required" << std::endl;
    consensus_print_help();
    return 1;
  }
  if (opt.bam.empty()) {
    std::cerr << "-b/--bam is required" << std::endl;
    consensus_print_help();
    return 1;
  }
  if (opt.allow_nonoverlapping_pair) {
    if (opt.pair_min_overlap != 0) {
      std::cerr << "-p/--pair_min_overlap has to be 0 if -i/--allow_nonoverlapping_pair is true\n";
      return 1;
    }
  }

  typedef int TValue;
  typedef Score<TValue, ScoreMatrix<Dna5, Default> > TScoringScheme;
  int const gapOpenScore = -8;
  int const gapExtendScore = -8;
  TScoringScheme scoringScheme(gapExtendScore, gapOpenScore);

  //std::cout << "User defined matrix (also Dna5 scoring matrix)..." << std::endl;
  setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());
  //showScoringMatrix(scoringScheme);

  SeqLib::BWAWrapper  bwa;
  bwa.LoadIndex(opt.reference);

  int8_t other_mm = -4;
  int8_t BS_mm = 0;
  int8_t mat[25] = {
    1, other_mm, BS_mm, other_mm, other_mm,
            other_mm, 1, other_mm, BS_mm, other_mm,
            other_mm, other_mm, 1, other_mm, other_mm,
            other_mm, other_mm, other_mm, 1, other_mm,
            other_mm, other_mm, other_mm, other_mm, other_mm
  };

  SeqLib::BamWriter bam_writer;
  SeqLib::BamWriter bam_e_writer;
  SeqLib::BamWriter bam_c_writer;
  string hdr_line;
  for (unsigned i = 0; i < bwa.GetIndex()->bns->n_seqs; ++i) {
    char hdr_line1[80];
    sprintf(hdr_line1, "@SQ\tSN:%s\tLN:%d", bwa.GetIndex()->bns->anns[i].name, bwa.GetIndex()->bns->anns[i].len);
    char hdr_line2[7];
    if (bwa.GetIndex()->bns->anns[i].is_alt) sprintf(hdr_line2, "\tAH:*\n");
    else sprintf(hdr_line2, "\n");
    char* full_text;
    full_text= (char*) malloc(strlen(hdr_line1)+strlen(hdr_line2)+1);
    strcpy(full_text, hdr_line1);
    strcat(full_text, hdr_line2);
    hdr_line += std::string(full_text);
    free(full_text);
  }

  bam_writer.Open(opt.outprefix + ".paired_reads.bam");
  bam_e_writer.Open(opt.outprefix + ".extended_reads.bam");
  bam_c_writer.Open(opt.outprefix + ".converted_reads.bam");
  std::vector<std::string> ref_dict;
  SeqLib::BamHeader bh(hdr_line);
  bam_writer.SetHeader(bh);
  bam_writer.WriteHeader();
  bam_e_writer.SetHeader(bh);
  bam_e_writer.WriteHeader();
  bam_c_writer.SetHeader(bh);
  bam_c_writer.WriteHeader();

  //std::cout << hdr_line << std::endl;
  cpputil::InsertSeqFactory isf(opt.bam,
                                opt.mapq,
                                opt.load_supplementary,
                                opt.load_secondary,
                                true,
                                !opt.allow_nonoverlapping_pair,
                                opt.clip3);
  int64_t read_counter = 0;
  int64_t good_read_counter = 0;
    while (!isf.finished()) {
      std::vector<cpputil::Segments> frag;
      while( (frag= isf.FetchReadNameSorted(true, false)).size() > 0) {
        for (auto seg : frag) {
          assert (seg.size() == 2);
          std::string read2 = seg[1].Sequence();
          cpputil::reverse_complement(read2);
//          int ol = cpputil::GetNumOverlapBasesPEAlignment(seg);
//          if (ol < opt.pair_min_overlap) continue;
          TSequence ref = seg[0].Sequence();
          TSequence query = read2;
          std::string converted_strand;
          std::string extended_strand;
          TAlign align1;
          seqan::resize(rows(align1), 2);
          seqan::assignSource(row(align1, 0), ref);
          seqan::assignSource(row(align1, 1), query);
          int s1 = seqan::globalAlignment(align1, scoringScheme,  seqan::AlignConfig<true, false, true, false>(), seqan::AffineGaps());
          if (s1 > opt.min_overlap) {
            ++good_read_counter;
            int a,b,c,d;
            std::tie(a,b,c,d) = print_AG_CT_mismatch(align1);
            float pi = 1.0 - (float) c/ (a+b+c+d);
            //std::cout << "l1: " <<seqan::length(ref) << " l2: " <<seqan::length(query) <<  " pi: " << pi <<" C/T: " << a <<" A/G " << b << std::endl;
            //std::cout << align1;
            bool first_read_extended;
            uint8_t *ext_qual = 0;
            uint8_t *cvt_qual = 0;
            if (a > b ) {
              //writer.WriteRecord(seg[0]);
              first_read_extended = true;
              extended_strand = seg[0].Sequence();
              converted_strand = seg[1].Sequence();
              ext_qual = bam_get_qual(seg[0].raw());
              cvt_qual = bam_get_qual(seg[1].raw());
            } else {
              first_read_extended = false;
              //writer.WriteRecord(seg[1]);
              extended_strand = seg[1].Sequence();
              converted_strand = seg[0].Sequence();
              ext_qual = bam_get_qual(seg[1].raw());
              cvt_qual = bam_get_qual(seg[0].raw());
            }
            mem_alnreg_v ar_et;
            ar_et = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                            extended_strand.length(),extended_strand.data());
            if (ar_et.n == 0 ) {
              std::cout << seg[0].Qname() << " has no alignment\n";
              continue;
            }
            size_t f_pidx = 0;
            for (size_t idx = 0; idx < ar_et.n; ++idx) {
              if (ar_et.a[idx].secondary < 0) {
                f_pidx = idx;
                break;
              }
            }

            mem_aln_t et_aln;

            int64_t rb=0, re=0;
            int rid;
            et_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, extended_strand.length(), extended_strand.c_str(), &ar_et.a[f_pidx]);
            rb = ar_et.a[f_pidx].rb - 200;
            re = ar_et.a[f_pidx].re + 200;
            int64_t l_pac = bwa.GetIndex()->bns->l_pac;
            if (re > l_pac <<1) re = l_pac<<1;
//            std::cout << "alignment on extended strand\n";
//            std::cout << bam_et;

            uint8_t *rseq = 0;
            rseq = bns_fetch_seq(bwa.GetIndex()->bns, bwa.GetIndex()->pac, &rb, (rb+re)>>1, &re, &rid);
            assert(rid == ar_et.a[f_pidx].rid);
            kswr_t ksw_aln;
            mem_alnreg_t cs_b;
            cpputil::reverse_complement(converted_strand);
            int l_ms = converted_strand.length();
            uint8_t* seq = 0;
            seq = (uint8_t *) malloc(l_ms);
            for (int i = 0; i < l_ms; ++i) // on the forward strand
              seq[i] = nst_nt4_table[(int)converted_strand[i]];
            int tmp, xtra = KSW_XSUBO | KSW_XSTART | (l_ms* bwa.GetMemOpt()->a < 250? KSW_XBYTE : 0) | (bwa.GetMemOpt()->min_seed_len * bwa.GetMemOpt()->a);

//            printf("* Global query:   "); for (int i = 0; i < l_ms; ++i) putchar("ACGTN"[(int)seq[i]]); putchar('\n');
//            printf("* Global ref: "); for (int i = 0; i < re-rb; ++i) putchar("ACGTN"[(int)rseq[i]]); putchar('\n');
            ksw_aln = ksw_align2(l_ms, seq, re - rb, rseq, 5, mat, bwa.GetMemOpt()->o_del, bwa.GetMemOpt()->e_del, bwa.GetMemOpt()->o_ins, bwa.GetMemOpt()->e_ins, xtra, 0);
            free(rseq);
            free(seq);
//            std::cout <<"score: "<< ksw_aln.score << std::endl;
            if (ksw_aln.score > opt.min_overlap) {
              cpputil::reverse_complement(converted_strand);
              cs_b.rid = rid;
              cs_b.is_alt = ar_et.a[f_pidx].is_alt;
              cs_b.qb = l_ms - (ksw_aln.qe + 1);
              cs_b.qe = l_ms - ksw_aln.qb;
              cs_b.rb = (l_pac<<1) - (rb + ksw_aln.te + 1);
              cs_b.re = (l_pac<<1) - (rb + ksw_aln.tb);
              cs_b.score = ksw_aln.score;
              cs_b.csub = ksw_aln.score2;
              cs_b.secondary = -1;
              cs_b.seedcov = (cs_b.re - cs_b.rb < cs_b.qe - cs_b.qb? cs_b.re - cs_b.rb : cs_b.qe - cs_b.qb) >> 1;
              auto cs_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, converted_strand.length(), converted_strand.c_str(), &cs_b);
              et_aln.flag |= 1;
              cs_aln.flag |= 1;

              if ((et_aln.is_rev) ^ (cs_aln.is_rev)) {
                et_aln.flag |= 2;
                cs_aln.flag |= 2;
              }
              if (cs_aln.is_rev) {
                et_aln.flag |= 0x20;
              }
              if (et_aln.is_rev) {
                cs_aln.flag |= 0x20;
              }
              if (first_read_extended) {
                et_aln.flag |= 0x40;
                cs_aln.flag |= 0x80;
              } else {
                et_aln.flag |= 0x80;
                cs_aln.flag |= 0x40;
              }
//              std::cout << "alignment on converted strand\n";
//              std::cout << bam_cs;
              auto bam_et = cpputil::BwaAlignment2BamRecord(et_aln, seg[0].Qname(), extended_strand, ext_qual);
              auto bam_cs = cpputil::BwaAlignment2BamRecord(cs_aln, seg[0].Qname(), converted_strand, cvt_qual);
              if (bam_cs.GetCigar().NumQueryConsumed() != converted_strand.length()) {
                std::cout << seg[0].Qname() << " cigar not interpretable\n";
                continue;
              }
              

              //add other but necessary infomation
              bam_et.shared_pointer()->core.mtid = bam_cs.ChrID();
              bam_et.shared_pointer()->core.mpos = bam_cs.Position();
              bam_cs.shared_pointer()->core.mtid = bam_et.ChrID();
              bam_cs.shared_pointer()->core.mpos = bam_et.Position();
              if (bam_et.Position() < bam_cs.Position()) {
                bam_et.shared_pointer()->core.isize = bam_cs.PositionEnd() - bam_et.Position() + 1;
                bam_cs.shared_pointer()->core.isize = - bam_cs.PositionEnd() + bam_et.Position() - 1;
              } else {
                bam_et.shared_pointer()->core.isize = -bam_et.PositionEnd() + bam_cs.Position() - 1;
                bam_cs.shared_pointer()->core.isize = bam_et.PositionEnd() - bam_cs.Position() + 1;
              }
              bam_writer.WriteRecord(bam_et);
              bam_e_writer.WriteRecord(bam_et);
              bam_writer.WriteRecord(bam_cs);
              bam_c_writer.WriteRecord(bam_cs);
              free(cs_aln.cigar);
            } else {
//              if (first_read_extended) et_aln.flag |= 1 | 0x40 | 8;
//              else et_aln.flag |= 1 | 0x80 | 8;
//              auto bam_et = cpputil::BwaAlignment2BamRecord(et_aln, seg[0].Qname(), extended_strand, ext_qual);
//              bam_writer.WriteRecord(bam_et);
            }
            free(ar_et.a);
            free(et_aln.cigar);
          }
          ++ read_counter;
        }
        if (read_counter == opt.max_read) {
          break;
        }
      }
      if (read_counter == opt.max_read) {
        break;
      }
    }
  std::cout << "generate " << read_counter << " consensus reads, where "<<good_read_counter <<" are good" << std::endl;
  return 0;
}
