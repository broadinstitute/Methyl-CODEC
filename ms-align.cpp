//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "spoa/spoa.hpp"

#include "BamIO.h"
#include "seqan/align.h"
#include "seqan_init_score.h"
#include "SeqLib/BWAWrapper.h"
#include "BamRecordExt.h"
#include "FastxRecord.h"
#include "FastxIO.h"
#include "AlignmentExtAndMethyl.h"
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
#include "DNAUtils.h"

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
  string itermol_out;
  int mapq = 10;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool clip3 = false;
//  int consensus_mode = 0;
  int pair_min_overlap = 1;
  int max_read = 0;
  int min_overlap = 30;
  bool call_overhang = false;
  string reference = "";
  string read_group_header = "";
//  bool strand_specific_fastq = false;
  int minbq = 0;
  int min_eof_dist = 0;
  string outprefix = "";
  int thread = 1;
  int MIN_READL = 15;
};

struct LibStat {
  int total_pairs = 0;
  int correct_paris = 0;
  int intermol_pairs = 0;
  int not_aligned_or_other_issues = 0;
};

static struct option  consensus_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"reference",                required_argument,      0,        'r'},
    {"load_supplementary",       no_argument,            0,        'l'},
    {"clip3",                    no_argument,            0,        'C'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"baseq",                    required_argument ,     0,        'q'},
    {"eof",                      required_argument ,     0,        'd'},
    {"outprefix",                required_argument ,     0,        'o'},
    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"call_overhang",            no_argument,            0,        't'},
//    {"strand_specific_fastq",   no_argument,            0,        's'},
    {"intermol_bam",            required_argument,            0,        'i'},
    {"max_read",                 required_argument,      0,        'M'},
    {"read_group_header",        required_argument,      0,        'R'},
    {0,0,0,0}
};

const char* consensus_short_options = "b:m:M:o:lCp:q:d:ti:r:R:";

void consensus_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              Input bam [required]\n";
  std::cerr<< "-o/--outprefix,                        Output sample prefix [required].\n";
  std::cerr<< "-m/--mapq,                             Min mapping quality [10].\n";
  std::cerr<< "-q/--baseq,                            Min base quality for calling both strand for calling metC[0].\n";
  std::cerr<< "-e/--min_eof_dist,                     Min distance to the end of the fragments for calling metC. [0].\n";
  std::cerr<< "-r/--reference,                        Reference for alignment. [null].\n";
  //Supplementary not support currently
  //std::cerr<< "-l/--load_supplementary,               Include supplementary alignment [false].\n";
//  std::cerr<< "-t/--call_overhang,                    Call metC in the overhang regions [false].\n";
  std::cerr<< "-C/--clip3,                            trim the 3'end soft clipping [false].\n";
//  std::cerr<< "-s/--strand_specific_fastq,            output fastq for each strand [false].\n";
//  std::cerr<< "-s/--output_singleend,                 The R1R2 consensus will be output in a single end format [false].\n";
  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair [1].\n";
//  std::cerr<< "-d/--dirtmp,                           Temporary dir for sorted bam [/tmp]\n";
//  std::cerr<< "-T/--thread,                           Number of threads for sort [1]\n";
  std::cerr<< "-i/--itermol_out,                      Output of non overlapping read-pairs, usually caused by intermolecular ligation. Bam output for protected strand and Fastq for original strand.  [False]\n";
  std::cerr<< "-R/--read_group_header,                read group header line such as '@RG\tID:foo\tSM:bar'. just like bwa\n";
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
      case 'R':
        opt.read_group_header = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.minbq = atoi(optarg);
        break;
      case 'd':
        opt.min_eof_dist = atoi(optarg);
        break;
      case 'l':
        opt.load_supplementary = true;
        break;
      case 't':
        opt.call_overhang = true;
        break;
      case 'i':
        opt.itermol_out = optarg;
        break;
      case 'C':
        opt.clip3 = true;
        break;
//      case 's':
//        opt.strand_specific_fastq = true;
//        break;
      case 'M':
        opt.max_read = atoi(optarg);
        break;
      case 'p':
        opt.pair_min_overlap = atoi(optarg);
        break;
      default:consensus_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


int codec_ms_align(int argc, char ** argv) {
  CssOptions opt;
  LibStat libstat;
  int parse_ret =  consensus_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    consensus_print_help();
    return 1;
  }
 // printf("%s\n", opt.read_group_header.c_str());
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

  int const gapOpenScore = -8;
  int const gapExtendScore = -8;
  cpputil::TScoringScheme scoringScheme(gapExtendScore, gapOpenScore);

  //std::cout << "User defined matrix (also Dna5 scoring matrix)..." << std::endl;
  setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());
  //showScoringMatrix(scoringScheme);

  SeqLib::BWAWrapper  bwa;
  bwa.LoadIndex(opt.reference);

  int8_t other_mm = -4;
  int8_t BS_mm = 0;
  int8_t N_mm = -1;
  int8_t mat[25] = {
    1, other_mm, BS_mm, other_mm, N_mm,
            other_mm, 1, other_mm, BS_mm, N_mm,
            other_mm, other_mm, 1, other_mm, N_mm,
            other_mm, other_mm, other_mm, 1, N_mm,
            N_mm, N_mm, N_mm, N_mm, N_mm
  };

  SeqLib::BamWriter bam_writer;
  cpputil::FastqWriter single_orig_strand_fq_writer;
  cpputil::FastqWriter single_prot_strand_fq_writer;
//  SeqLib::BamWriter bam_e_writer;
//  SeqLib::BamWriter bam_c_writer;
//  cpputil::FastqWriter fq_e_writer, fq_c_writer;
//  if (opt.strand_specific_fastq) {
//    fq_e_writer.open(opt.outprefix + ".protected_reads.fastq.gz");
//    fq_c_writer.open(opt.outprefix + ".converted_reads.fastq.gz");
//  }
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
  if (!opt.read_group_header.empty()) {
    size_t pos = 0;
    while ((pos = opt.read_group_header.find("\\t", pos)) != std::string::npos) {
      opt.read_group_header.replace(pos, 2, "\t");
      pos += 1; // Move past the replaced character
    }
    hdr_line += opt.read_group_header + "\n";
  }

  bam_writer.Open(opt.outprefix + ".paired_reads.bam");

  std::vector<std::string> ref_dict;
  SeqLib::BamHeader bh(hdr_line);
  bam_writer.SetHeader(bh);
  bam_writer.WriteHeader();

  if (not opt.itermol_out.empty()) {
    single_orig_strand_fq_writer.open(opt.itermol_out + ".converted_strand.fastq.gz");
    single_prot_strand_fq_writer.open(opt.itermol_out + ".protected_strand.fastq.gz");
  }

  SeqLib::RefGenome refseq;
  refseq.LoadIndex(opt.reference);

  //std::cout << hdr_line << std::endl;
  cpputil::InsertSeqFactory isf(opt.bam,
                                opt.mapq,
                                opt.load_supplementary,
                                opt.load_secondary,
                                true,
                                false,
                                opt.clip3);
    while (!isf.finished()) {
      std::vector<cpputil::Segments> frag;
      while( (frag= isf.FetchReadNameSorted(true, false)).size() > 0) {
        for (auto seg : frag) {
          assert (seg.size() == 2);
          ++libstat.total_pairs;
          std::string rx1, rx2;
          if (not seg[0].GetZTag("RX", rx1)) throw std::runtime_error("No RX tag found\n");
          if (not seg[1].GetZTag("RX", rx2)) throw std::runtime_error("No RX tag found\n");
          int is_correct_ms =cpputil::IsCorrectMSPairedReads2(seg, scoringScheme, opt.pair_min_overlap);
          if (is_correct_ms > 0) {
            std::string converted_strand;
            std::string extended_strand;
            bool first_read_extended;
            uint8_t *ext_qual = 0;
            uint8_t *cvt_qual = 0;
            if (is_correct_ms == 1) {
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

//          alignment on extended strand using BWA-mem
            mem_alnreg_v ar_et;
            ar_et = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                            extended_strand.length(),extended_strand.data());
            if (ar_et.n == 0 ) {
              //std::cout << seg[0].Qname() << " has no alignment\n";
              ++libstat.not_aligned_or_other_issues;
              continue;
            }
            size_t f_pidx = 0, f_2idx;
            int sec_as = 0;
            for (size_t idx = 0; idx < ar_et.n; ++idx) {
              if (ar_et.a[idx].secondary < 0) {
                f_pidx = idx;
              } else {
                sec_as = std::max(sec_as, ar_et.a[idx].score);
              }
            }
            mem_aln_t et_aln;
            et_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, extended_strand.length(), extended_strand.c_str(), &ar_et.a[f_pidx]);

//          map converted strand using BWA-SW at position where extended strand was mapped
            int64_t rb=0, re=0;
            rb = ar_et.a[f_pidx].rb - 500;
            re = ar_et.a[f_pidx].re + 500;
            int64_t l_pac = bwa.GetIndex()->bns->l_pac;
            if (re > l_pac <<1) re = l_pac<<1;
            uint8_t *rseq = 0;
            int rid;
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
            if (ksw_aln.score > opt.min_overlap) { // if protected strand guided alignment is successful
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
              bam_et.AddIntTag("XC", 0);
              auto bam_cs = cpputil::BwaAlignment2BamRecord(cs_aln, seg[0].Qname(), converted_strand, cvt_qual);
              bam_cs.AddIntTag("XC", 1);
              bam_et.AddIntTag("MQ", bam_cs.MapQuality());
              bam_cs.AddIntTag("MQ", bam_et.MapQuality());
              bam_et.AddZTag("MC", bam_cs.CigarString());
              bam_cs.AddZTag("MC", bam_et.CigarString());
              if (bam_cs.GetCigar().NumQueryConsumed() != converted_strand.length()) {
                std::cerr << seg[0].Qname() << " cigar not interpretable\n";
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
              cpputil::Segments aligned_segs = {bam_et, bam_cs};
              auto metc = cpputil::CallingMetC(refseq, bh, aligned_segs, opt.call_overhang, opt.minbq, opt.min_eof_dist);

              bam_et.AddIntTag("XS", sec_as);
              bam_cs.AddZTag("XM", metc);
              bam_cs.AddZTag("XR", "GA");
              if (bam_cs.ReverseFlag()) {
                bam_cs.AddZTag("XG", "CT");
              } else {
                bam_cs.AddZTag("XG", "GA");
              }
              if (first_read_extended) {
                bam_et.AddZTag("RX", rx1);
                bam_cs.AddZTag("RX", rx2);
              } else {
                bam_et.AddZTag("RX", rx2);
                bam_cs.AddZTag("RX", rx1);
              }
              bam_writer.WriteRecord(bam_et);
              bam_writer.WriteRecord(bam_cs);

              //make a single-end bam record of converted strand for bismark_extract_methylation
//              if (bam_cs.ReverseFlag()) {
//                bam_cs.shared_pointer()->core.flag = 16;
//              } else {
//                bam_cs.shared_pointer()->core.flag = 0;
//              }
//              bam_cs.shared_pointer()->core.mtid = -1;
//              bam_cs.shared_pointer()->core.mpos = -1;
//              bam_cs.shared_pointer()->core.isize = 0;
//              metbam_writer.WriteRecord(bam_cs);

              //output fastq for each strands
//              if (opt.strand_specific_fastq) {
//                cpputil::FastxRecord fq_et(bam_et);
//                fq_e_writer.Write(fq_et);
//                cpputil::FastxRecord fq_cs(bam_cs);
//                fq_c_writer.Write(fq_cs);
//              }
              ++libstat.correct_paris;
              free(cs_aln.cigar);
            } else { // converted strand failed to align to the same place where extended strand got aligned
              if (not opt.itermol_out.empty()) {
                std::string qname = seg[0].Qname();
                if (first_read_extended) {
                  single_prot_strand_fq_writer.Write(qname +"/1", seg[0].Sequence(), seg[0].Qualities());
                  single_orig_strand_fq_writer.Write(qname + "/2", seg[1].Sequence(), seg[1].Qualities());
                } else {
                  single_prot_strand_fq_writer.Write(qname + "/2", seg[1].Sequence(), seg[1].Qualities());
                  single_orig_strand_fq_writer.Write(qname + "/1", seg[0].Sequence(), seg[0].Qualities());
                }
              }
              ++libstat.intermol_pairs;
            }
            free(ar_et.a);
            free(et_aln.cigar);
          } else { // not correct MS reads
            std::string r1 = seg[0].Sequence();
            std::string r2 = seg[1].Sequence();
            if (r1.length() > opt.MIN_READL and r2.length() > opt.MIN_READL) {
              ++libstat.intermol_pairs;
            }
            if (not opt.itermol_out.empty()) {
              std::string q1 = seg[0].Qualities();
              std::string q2 = seg[1].Qualities();
              std::string name = seg[0].Qname();

              if (cpputil::is_bisulfite_converted(r1, opt.MIN_READL)) {
                //Need to replace with bamwriter so we can keep the RX tag
                single_orig_strand_fq_writer.Write(name + "/1", r1, q1);
              } else {
                single_prot_strand_fq_writer.Write(name +"/1", r1, q1);
              }

              if (cpputil::is_bisulfite_converted(r2, opt.MIN_READL)) {
                //Need to replace with bamwriter so we can keep the RX tag
                single_orig_strand_fq_writer.Write(name + "/2", r2, q2);
              } else {
                single_prot_strand_fq_writer.Write(name + "/2", r2, q2);
              }
            }
          }
        }
        if (libstat.total_pairs == opt.max_read) {
          break;
        }
      }
      if (libstat.total_pairs == opt.max_read) {
        break;
      }
    }
  std::cout << "total\tcorrect\tintermolecular\tnot_aligned_or_other_issues\n";
  std::cout << libstat.total_pairs <<"\t" << libstat.correct_paris << "\t" << libstat.intermol_pairs << "\t" << libstat.not_aligned_or_other_issues << std::endl;
  return 0;
}
