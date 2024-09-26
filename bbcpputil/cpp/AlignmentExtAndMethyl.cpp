//
// Created by Ruolin Liu on 3/19/20.
//

#include "AlignmentExtAndMethyl.h"
#include "DNAUtils.h"

namespace cpputil {

std::string GetConsensusTemplate(const Segments& segs, int32_t& ref_most_left) {
  /*
   *  '.' : uninitialized
   *  '+' : insertion
   */
  ref_most_left = std::numeric_limits<int32_t>::max();
  int32_t ref_most_right = 0;
  for (auto & seg : segs) {
    ref_most_left = std::min(seg.PositionWithSClips(), ref_most_left);
    ref_most_right = std::max(seg.PositionEndWithSClips(), ref_most_right);
  }
  int32_t ref_span = ref_most_right - ref_most_left;
  assert(ref_span > 0);
  std::map<int,int> ins_len;
  for (auto & seg: segs) {
    const SeqLib::Cigar cigar = seg.GetCigar();
    int32_t relative_start = seg.PositionWithSClips() - ref_most_left;
    find_insert_(cigar, relative_start, ins_len);
  }
//  const SeqLib::Cigar left_cigar = seg.front().GetCigar();
//  const SeqLib::Cigar right_cigar = seg.back().GetCigar();
//  find_insert_(left_cigar, 0, ins_len);
//  find_insert_(right_cigar, seg.back().PositionWithSClips() - seg.front().PositionWithSClips(), ins_len);
  std::string consens;
  int b = 0;
  for (const auto pos_len : ins_len) {
    consens += std::string(pos_len.first - b, '.');
    consens += std::string(pos_len.second, '+');
    b = pos_len.first;
  }
  consens += std::string(ref_span - b, '.');
  return consens;
}

std::string GetGappedSeqForRef(const SeqLib::RefGenome&ref, const SeqLib::BamHeader& header, const SeqLib::BamRecord &r, const int start, const std::string& consensus_template) {
  // seq should be same length as r.Sequence() but can with difference sequence
  // start should be relative start position to the fragment
  // `.` is overhang
  // `+` is insertion
  // `-` is deletion
  int template_start = 0;
  std::string consns_templ = consensus_template;
  auto seq = ref.QueryRegion(header.IDtoName(r.ChrID()), r.PositionWithSClips(), r.PositionEndWithSClips());
  int ref_count = 0;
  for (char c : consns_templ) {
    if (ref_count == start) break;
    template_start++;
    if (c == '.') ref_count++;
  }
  auto const &cigar = r.GetCigar();

  int ref_start = 0;
  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    if (c->Type() == 'S' || c->Type() == 'M' || c->Type() == 'D') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        while ('+' == consns_templ[template_start]) { template_start++; };
        consns_templ[template_start++] = toupper(seq[ref_start++]);
      }
    } else if (c->Type() == 'I') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        consns_templ[template_start++] = '+';
      }
      while (consns_templ[template_start] == '+') { template_start++; };
    }
  }
  return consns_templ;
}

std::pair<std::string, std::string>
    GetGappedSeqAndQual(const SeqLib::BamRecord &r, const std::string& seq, const int start, const std::string& consensus_template, bool include_soft_clip = true ) {

  // seq should be same length as r.Sequence() but can with difference sequence
  // `.` is overhang
  // `+` is insertion
  // `-` is deletion
  int template_start = 0;
  std::string consns_templ = consensus_template;
  std::string quality_templ = std::string(consensus_template.size(), 33);
  int ref_count = 0;
  for (char c : consns_templ) {
    if (ref_count == start) break;
    template_start++;
    if (c == '.') ref_count++;
  }
  auto const &cigar = r.GetCigar();
  auto const &qual = r.Qualities();

  int read_start = 0;
  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    if (c->Type() == 'S' || c->Type() == 'M' || c->Type() == 'D') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        while ('+' == consns_templ[template_start]) { template_start++; };
        char s, q;
        if (c->Type() == 'D') {
          s = '-';
          q = qual[read_start];
        } else if (c->Type () == 'S' and not include_soft_clip) {
          s = 'N';
          q = qual[read_start++];
        } else {
          s = seq[read_start];
          q = qual[read_start++];
        }
        if (consns_templ[template_start] != '+') {
          consns_templ[template_start] = s;
          quality_templ[template_start++] = q;
        }
      }
    } else if (c->Type() == 'I') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        consns_templ[template_start] = seq[read_start];
        quality_templ[template_start++] = qual[read_start++];
      }
      while (consns_templ[template_start] == '+') { template_start++; };
    }
  }
  return std::make_pair(consns_templ, quality_templ);
}

std::string MergePairSeq(const Segments &seg, const std::vector<std::string>& seqs, bool trim_overhang) {
  //seqs should hold the fastq seq for seg. They could be same length but different seqs
  //not to change indel
  assert (seg.size() == 2);
  int ref_most_left;
  const char NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(seg, ref_most_left);
  std::string consns_seq1(consns_templ.size(), NUL);
  int32_t start = 0;
  std::string dnaseq, qual;
  std::vector<std::string> dna_pileup(seg.size());
  std::vector<std::string> qual_pileup(seg.size());
  for (unsigned sid = 0; sid < seg.size(); ++sid) {
    auto const& piece = seg[sid];
    start = piece.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(piece, seqs[sid], start, consns_templ);
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
  for (unsigned jj = 0; jj < consns_templ.size(); ++jj) {
    // paired baseq calibration. If only one of the baseq < cutoff, make the other one baseq = cutoff -1
    // so that when later we filter by baseq by this cutoff, they either both stay or both out
    if (dna_pileup[0][jj] == '.' or dna_pileup[1][jj] == '.') { // overhang
      if (not trim_overhang) {
        if (dna_pileup[0][jj] != '.' and dna_pileup[0][jj] != '-' and dna_pileup[1][jj] == '.') {
          consns_seq1[jj] = dna_pileup[0][jj];
        } else if (dna_pileup[0][jj] == '.' and dna_pileup[1][jj] != '.' and dna_pileup[1][jj] != '-') {
          consns_seq1[jj] = dna_pileup[1][jj];
        }
      }
    } else {
      if (dna_pileup[0][jj] != dna_pileup[1][jj]) {
        assert(dna_pileup[0][jj] != '-' or dna_pileup[1][jj] != '+');
        assert(dna_pileup[0][jj] != '+' or dna_pileup[1][jj] != '-');
        if (dna_pileup[0][jj] == '-' or dna_pileup[0][jj] == '+') {
          consns_seq1[jj] = dna_pileup[1][jj];
        }
        else if (dna_pileup[1][jj] == '-' or dna_pileup[1][jj] == '+') {
          consns_seq1[jj] = dna_pileup[0][jj];
        }
        else {
          if (qual_pileup[0][jj] < qual_pileup[1][jj]) {
            consns_seq1[jj] = dna_pileup[1][jj];
          } else {
            consns_seq1[jj] = dna_pileup[0][jj];
          }
        }

      } else if (dna_pileup[0][jj] >= 'A') {
        consns_seq1[jj] = dna_pileup[0][jj];
      } else if (dna_pileup[0][jj] == '+') {
        assert(false);
      }
    }
  }
  consns_seq1.erase(std::remove(consns_seq1.begin(), consns_seq1.end(), NUL), consns_seq1.end());
  return consns_seq1;
}

std::string MergePair(const Segments &seg, bool trim_overhang) {
  std::vector<std::string> seqs;
  std::vector<std::string> dummy_quals;
  for (auto&s : seg) {
    seqs.push_back(s.Sequence());
  }
  auto seq = MergePairSeq(seg, seqs, trim_overhang);
  return seq;
}

//std::pair<std::string, std::string> PairSeqConsensus(const Segments &seg, bool trim_overhang, int qcutoff) {
//  std::vector<std::string> seqs;
//  std::vector<std::string> dummy_quals;
//  for (auto&s : seg) {
//    seqs.push_back(s.Sequence());
//  }
//  auto seq = cpputil::PairConsensus(seg, seqs, trim_overhang, qcutoff, dummy_quals);
//  return seq;
//}

std::pair<std::vector<std::string>, std::vector<std::string>> GetPairPileup(const Segments &seg) {
  assert(seg.size() == 2);
  const char NUL = 6;
  int ref_most_left;
  const std::string consns_templ = GetConsensusTemplate(seg, ref_most_left);
  std::vector<std::string> seqs;
  for (auto&s : seg) {
    seqs.push_back(s.Sequence());
  }
  std::string consns_seq1(consns_templ.size(), NUL);
  std::string consns_seq2(consns_templ.size(), NUL);
  int32_t start = 0;
  std::string dnaseq, qual;
  std::vector<std::string> dna_pileup(seg.size());
  std::vector<std::string> qual_pileup(seg.size());
  for (unsigned sid = 0; sid < seg.size(); ++sid) {
    auto const& piece = seg[sid];
    start = piece.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(piece, seqs[sid], start, consns_templ);
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
  return std::make_pair(dna_pileup, qual_pileup);
}

char ResolveCytosineContextSS(const std::string& cvt_dna_str, const std::string& qual_str, const std::string& ref_str, int jj, int qcutoff_char, bool cvs_is_reverse) {
    if ((cvt_dna_str[jj] == ref_str[jj])  or
        (cvs_is_reverse && cvt_dna_str[jj] == 'T' and ref_str[jj] == 'C') or
    (not cvs_is_reverse && cvt_dna_str[jj] == 'A' and ref_str[jj] == 'G'))
    {
      return ref_str[jj];
    }
    if (qual_str[jj] >= qcutoff_char) {
      return cvt_dna_str[jj];
    } else {
      return ref_str[jj];
    }
}

char ResolveCytosineContextDS(const std::string& cvt_dna_str, const std::string& cvt_qual_str,
                              const std::string& prt_dna_str, const std::string& prt_qualt_str,
                              const std::string& ref_str, int jj, int qcutoff_char, bool cvs_is_reverse) {
  if (cvt_qual_str[jj] >= qcutoff_char and prt_qualt_str[jj] >= qcutoff_char) {
    if (cvt_dna_str[jj] == prt_dna_str[jj]) {
      return prt_dna_str[jj];
    } else if ((cvs_is_reverse and prt_dna_str[jj] == 'C' and cvt_dna_str[jj] == 'T') or
            (not cvs_is_reverse and prt_dna_str[jj] == 'G' and cvt_dna_str[jj] == 'A')) {
      return prt_dna_str[jj];
    } else {
      return 'N';
    }
  } else {
    if (cvt_dna_str[jj] == prt_dna_str[jj] and prt_dna_str[jj] == ref_str[jj]) {
      return ref_str[jj];
    } else if ((cvs_is_reverse and prt_dna_str[jj] == 'C' and cvt_dna_str[jj] == 'T' and ref_str[jj] == 'C') or
               (not cvs_is_reverse and prt_dna_str[jj] == 'G' and cvt_dna_str[jj] == 'A' and ref_str[jj] == 'G')) {
      return ref_str[jj];
    } else {
      return 'N';
    }
  }
}

std::string CallMetC(const SeqLib::RefGenome& ref, const SeqLib::BamHeader& bamheader, const Segments &segs, bool call_overhang, int qcutoff, int eof) {
  // eof: min distance to the end of fragment filter
  // this will not call methylation status at soft clip regions by setting include_soft_clip = false in GetGappedSeqAndQual() which return 'N' at soft clip regions
  std::string metstr;
  assert (segs.size() == 2);
  int32_t ref_most_left;
  const char NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  if (ref_most_left < 0) {
    return metstr;
  }
  int32_t start = 0;
  std::string dnaseq, qual, ref_pileup;
  std::vector<std::string> dna_pileup(segs.size());
  std::vector<std::string> qual_pileup(segs.size());
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seg.Sequence(), start, consns_templ, false);
    if (sid == 1) { // ref pileup only for converted strand
      ref_pileup = GetGappedSeqForRef(ref, bamheader, seg, start, consns_templ);
    }
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
//  int r1_converted;
//  bool r1_getxc = segs[0].GetIntTag("XC", r1_converted);
//  int r2_converted;
//  bool r2_getxc = segs[1].GetIntTag("XC", r2_converted);
//  if (not r1_getxc | not r2_getxc | not (r1_converted ^ r2_converted)) {
//    std::cerr << "Methylation status unclear\n";
//    return metstr;
//  }
  int cs_idx = 1;
  int cs_reverse = segs[cs_idx].ReverseFlag();


  char qcutoff_char = static_cast<char>(33 + qcutoff);
  //unsigned refpos = 0;
  //char refc = 0;
  for (unsigned jj = 0; jj < consns_templ.size(); ++jj) {
    //if (consns_templ[jj] != '+') refc = toupper(refseq[refpos++]);
    if (jj < eof or jj + eof >= consns_templ.size()) {
      if (IsdNTP(dna_pileup[cs_idx][jj])) {
        metstr += ".";
      }
      continue;
    }
    if (dna_pileup[1-cs_idx][jj] == '.') { // overhang
      if (call_overhang) {
        if (cs_reverse) {
          if ((dna_pileup[cs_idx][jj] == 'C' or dna_pileup[cs_idx][jj] == 'T') and qual_pileup[cs_idx][jj] >= qcutoff_char) {
            CCTX cctx = CCTX::Other;
            if (jj + 1 < consns_templ.size()) {
              char c1 = ResolveCytosineContextSS(dna_pileup[cs_idx], qual_pileup[cs_idx], ref_pileup, jj+1, qcutoff_char, cs_reverse);
              if (c1 == 'G') {
                cctx = CCTX::CpG;
              } else if (IsH(c1)) {
                if (jj + 2 < consns_templ.size()) {
                  char c2 = ResolveCytosineContextSS(dna_pileup[cs_idx], qual_pileup[cs_idx], ref_pileup, jj+2, qcutoff_char, cs_reverse);
                  if (IsH(c2)) {
                    cctx = CCTX::CHH;
                  } else if (dna_pileup[cs_idx][jj+2] == 'G') {
                    cctx = CCTX::CHG;
                  }
                }
              }
            }

            switch (dna_pileup[cs_idx][jj]) {
              case 'C':
                metstr += std::string(1, BisMarkSymb(true, cctx));
                break;
              case 'T':
                metstr += std::string(1, BisMarkSymb(false, cctx));
                break;
            }
          } else if (dna_pileup[cs_idx][jj] >= 'A')  { // not a C
            metstr += ".";
          }
        } else { // on forward strand
          if ((dna_pileup[cs_idx][jj] == 'G' or dna_pileup[cs_idx][jj] == 'A') and qual_pileup[cs_idx][jj] >= qcutoff_char) {
            CCTX cctx = CCTX::Other;
            if (jj > 0) {
              char c1 = ResolveCytosineContextSS(dna_pileup[cs_idx], qual_pileup[cs_idx], ref_pileup, jj-1, qcutoff_char, cs_reverse);
              if (c1 == 'C') {
                cctx = CCTX::CpG;
              } else if (IsComplementH(c1)) {
                if (jj > 1) {
                  char c2 = ResolveCytosineContextSS(dna_pileup[cs_idx], qual_pileup[cs_idx], ref_pileup, jj-2, qcutoff_char, cs_reverse);
                  if (IsComplementH(c2)) {
                    cctx = CCTX::CHH;
                  } else if (c2 == 'C') {
                    cctx = CCTX::CHG;
                  }
                }
              }
            }

            switch (dna_pileup[cs_idx][jj]) {
              case 'G':
                metstr += std::string(1, BisMarkSymb(true, cctx));
                break;
              case 'A':
                metstr += std::string(1, BisMarkSymb(false, cctx));
                break;
            }
          } else if (dna_pileup[cs_idx][jj] >= 'A')  { // not a C
            metstr += ".";
          }
        }
      } else { // do not call overhang
        if (dna_pileup[cs_idx][jj] >= 'A') {
          metstr += ".";
        }
      }
    } else { // not overhang
      if (cs_reverse) { // cs reverse strand
        if (dna_pileup[1-cs_idx][jj] == 'C' and
            //ref_pileup[jj] == 'C' and
            qual_pileup[1-cs_idx][jj] >= qcutoff_char and
            qual_pileup[cs_idx][jj] >= qcutoff_char) {

          CCTX cctx = CCTX::Other;
          if (jj + 1 < dna_pileup[cs_idx].size()) {
            auto c1 = ResolveCytosineContextDS(dna_pileup[cs_idx], qual_pileup[cs_idx], dna_pileup[1-cs_idx], qual_pileup[1-cs_idx], ref_pileup, jj+1, qcutoff_char, cs_reverse);
            if (c1 == 'N') {
              cctx = CCTX::Other;
            }
            else if (c1 == 'G') {
              cctx = CCTX::CpG;
            } else if (IsH(c1)) {
              if (jj + 2 < dna_pileup[cs_idx].size()) {
                auto c2 = ResolveCytosineContextDS(dna_pileup[cs_idx], qual_pileup[cs_idx], dna_pileup[1-cs_idx], qual_pileup[1-cs_idx], ref_pileup, jj+2, qcutoff_char, cs_reverse);
                if (IsH(c2)) {
                  cctx = CCTX::CHH;
                } else if (c2 == 'G') {
                  cctx = CCTX::CHG;
                }
              }
            }
          }

          switch (dna_pileup[cs_idx][jj]) {
            case 'C':
              metstr += std::string(1, BisMarkSymb(true, cctx));
              break;
            case 'T':
              //if (refc == 'C') {
              metstr += std::string(1, BisMarkSymb(false, cctx));
//              } else {
//                metstr += ".";
//              }
              break;
            case 'G':
            case 'A': // seq error
            case 'N':
              metstr += ".";
              break;
            default:
              break;
          }
        } else if (dna_pileup[cs_idx][jj] >= 'A') {
          metstr += ".";
        }
      } else { // cs forward strand
        if (dna_pileup[1-cs_idx][jj] == 'G' and
            //ref_pileup[jj] == 'G' and
            qual_pileup[1-cs_idx][jj] >= qcutoff_char and
            qual_pileup[cs_idx][jj] >= qcutoff_char) {

          CCTX cctx = CCTX::Other;
          if (jj > 0) {
            auto c1 = ResolveCytosineContextDS(dna_pileup[cs_idx], qual_pileup[cs_idx], dna_pileup[1-cs_idx], qual_pileup[1-cs_idx], ref_pileup, jj-1, qcutoff_char, cs_reverse);
            if (c1 == 'N') {
              cctx = CCTX::Other;
            }
            else if (c1 == 'C') {
              cctx = CCTX::CpG;
            } else if (IsComplementH(c1)) {
              if (jj > 1) {
                auto c2 = ResolveCytosineContextDS(dna_pileup[cs_idx], qual_pileup[cs_idx], dna_pileup[1-cs_idx], qual_pileup[1-cs_idx], ref_pileup, jj-2, qcutoff_char, cs_reverse);
                if (c2 == 'N') {
                  cctx = CCTX::Other;
                }
                else if (IsComplementH(c2)) {
                  cctx = CCTX::CHH;
                } else if (c2 == 'C') {
                  cctx = CCTX::CHG;
                }
              }
            }
          }

          switch (dna_pileup[cs_idx][jj]) {
            case 'G':
              metstr += std::string(1, BisMarkSymb(true, cctx));
              break;
            case 'A':
              metstr += std::string(1, BisMarkSymb(false, cctx));
              break;
            case 'C':
            case 'T':
            case 'N':
              metstr += ".";
              break;
            default:
              break;
          }
        } else if (dna_pileup[cs_idx][jj] >= 'A') {
          metstr += ".";
        }
      }
    }
  }
  return metstr;
}


// Simple consensus by requiring all bases to be the same, otherwise 'N'.
// Same for INS and DEL. The Inserted seq has to be the same. The insertions with different lengths are truncated the the
// smallest length.
// The consensus of a deleted base and undeleted base is N
std::pair<std::string, std::string> MergeSegs(const Segments &segs, const std::vector<std::string> seqs,
                                              bool trim_overhang, int qcutoff,
                                              std::vector<std::string>& ori_quals) {
  //seqs should hold the fastq seq for segs. They could be same length but different seqs
  //Generate a consensus template by opening gaps(MSA format)
  //Each read in ungapped space will be converted to gap space and iteratively updates the consensus
  if (segs.empty()) {
    assert(false);
  }
  int ref_most_left;
  const int NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  std::string consns_seq(consns_templ.size(), NUL);
  std::string consns_qual(consns_templ.size(), 33);
  ori_quals.resize(segs.size());
  for (auto & qual : ori_quals) {
    qual = std::string(consns_templ.size(), 33);
  }

  int32_t start = 0;
  std::string dnaseq, qual;
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seqs[sid], start, consns_templ);
    for (unsigned ii = 0; ii < consns_templ.size(); ++ii) {
      if (trim_overhang && dnaseq[ii] == '.') {
        consns_seq[ii] = '+'; // to be removed
      }
      if (dnaseq[ii] == '.') continue;
      if (consns_seq[ii] == NUL) {
        consns_seq[ii] = dnaseq[ii];
      } else if (consns_seq[ii] != 'N' && consns_seq[ii] != '+') {
         if (dnaseq[ii] != consns_seq[ii]) {
           consns_seq[ii] = dnaseq[ii] == '+' ? '+' : 'N';
         }
      }
      consns_qual[ii] = std::max(qual[ii], consns_qual[ii]);
      ori_quals[sid][ii] = qual[ii];
    }
  }
  for (unsigned ii = 0; ii < consns_templ.size(); ++ii) {
    if (consns_seq[ii] == NUL || consns_seq[ii] == '+' || consns_seq[ii] == '-') {
      consns_seq[ii] = NUL;
      consns_qual[ii] = NUL;
      for (auto& qual : ori_quals) {
        qual[ii] = NUL;
      }
    }
  }
  consns_seq.erase(std::remove(consns_seq.begin(), consns_seq.end(), NUL), consns_seq.end());
  consns_qual.erase(std::remove(consns_qual.begin(), consns_qual.end(), NUL), consns_qual.end());
  for (auto& qual : ori_quals) {
    qual.erase(std::remove(qual.begin(),qual.end(), NUL ), qual.end());
  }
  return std::make_pair(consns_seq, consns_qual);
}

int IsCorrectMSPairedReads2(const Segments& seg, const TScoringScheme& ss, const int min_score) {
  /*
   * return 0: for incorrect MS reads
   * 1: correct MS reads with first read is the protected strand
   * 2: correct MS reads with second read is the protected strand
   */
  int ret = IsCorrectMSPairedReads1(seg);
  if (ret != 0) {
    return ret;
  }
  std::string read2 = seg[1].Sequence();
  cpputil::reverse_complement(read2);
  TSequence ref = seg[0].Sequence();
  TSequence query = read2;
  TAlign align1;
  seqan::resize(rows(align1), 2);
  seqan::assignSource(row(align1, 0), ref);
  seqan::assignSource(row(align1, 1), query);
  int s1 = seqan::globalAlignment(align1, ss,  seqan::AlignConfig<true, false, true, false>(), seqan::AffineGaps());
  if (s1 > min_score) {
    int a,b,c,d;
    std::tie(a,b,c,d) = print_AG_CT_mismatch(align1);
    return a>b ? 1 : 2;
    //float pi = 1.0 - (float) c/ (a+b+c+d);
    //std::cout << "l1: " <<seqan::length(ref) << " l2: " <<seqan::length(query) <<  " pi: " << pi <<" C/T: " << a <<" A/G " << b << std::endl;
    //std::cout << align1;

  } else {
    return 0;
  }
}

int GetConvertStrandIndex(const Segments& seg) {
  int xc1, xc2;
  int cidx = -1;
  bool xc1_flag = seg[0].GetIntTag("XC", xc1);
  bool xc2_flag = seg[1].GetIntTag("XC", xc2);
  if (xc1_flag and xc2_flag and xc1+xc2 == 1) {
    cidx = xc1 == 1? 0 : 1;
  }
  return cidx;
}

void ResolveCT_GA_bases_MSPairedReads(const SeqLib::RefGenome& ref, const SeqLib::BamHeader& header, std::vector<Segments>& frag) {
  for (auto & seg : frag) {
    int cidx = GetConvertStrandIndex(seg);
    if (cidx != -1) {
      int ref_most_left;
      const char NUL = 6;
      const std::string consns_templ = GetConsensusTemplate(seg, ref_most_left);
      //std::string cnvt_ref_templ = GetGappedSeqForRef(ref, header, seg[cidx], seg[cidx].PositionWithSClips() - ref_most_left, consns_templ);
      std::string cnvt_read_seq(consns_templ.size(), NUL);

      std::vector<std::string> dna_pileup, qual_pileup;
      std::tie(dna_pileup, qual_pileup) = GetPairPileup(seg);
      for (unsigned jj = 0; jj < dna_pileup[0].size(); ++jj) {
        if (dna_pileup[cidx][jj] >= 'A') {
           if (seg[cidx].ReverseFlag()) {
             if (dna_pileup[cidx][jj] == 'T' and dna_pileup[1-cidx][jj] == 'C' /*and cnvt_ref_templ[jj] != 'T'*/) {
               cnvt_read_seq[jj] = 'C';
             } else {
               cnvt_read_seq[jj] = dna_pileup[cidx][jj];
             }
           } else {
             if (dna_pileup[cidx][jj] == 'A' and dna_pileup[1-cidx][jj] == 'G' /*and cnvt_ref_templ[jj] != 'A'*/) {
               cnvt_read_seq[jj] = 'G';
             } else {
                cnvt_read_seq[jj] = dna_pileup[cidx][jj];
             }
           }
        }
      }
      cnvt_read_seq.erase(std::remove(cnvt_read_seq.begin(), cnvt_read_seq.end(), NUL), cnvt_read_seq.end());
      auto qual = seg[cidx].Qualities(0);
      seg[cidx].SetSequence(cnvt_read_seq);
      seg[cidx].SetQualities(qual, 0);
    }
  }
  return;
}

SeqLib::BamRecord SingleEndBWA(const SeqLib::BWAWrapper& bwa, const SeqLib::BamRecord& ubam, const int MIN_READL) {
  mem_alnreg_v ar1;
  mem_aln_t ar1_aln;
  int ar1_sec_as = 0;
  std::string r1 = ubam.Sequence();
  uint8_t *qual1 = 0;
  qual1 = bam_get_qual(ubam.raw());
  if (r1.length() < MIN_READL) {
    return ubam;
  } else {
    ar1 = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                     r1.length(),r1.data());
    if (ar1.n > 0 ) {
      size_t f_pidx = 0;
      for (size_t idx = 0; idx < ar1.n; ++idx) {
        if (ar1.a[idx].secondary < 0) {
          f_pidx = idx;
        } else {
          ar1_sec_as = std::max(ar1_sec_as, ar1.a[idx].score);
        }
      }
      ar1_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, r1.length(), r1.c_str(), &ar1.a[f_pidx]);
      ar1_aln.flag |= 1;
      if (ubam.FirstFlag()) {
        ar1_aln.flag |= 0x40;
      } else {
        ar1_aln.flag |= 0x80;
      }

      auto bam1 = cpputil::BwaAlignment2BamRecord(ar1_aln, ubam.Qname(), r1, qual1);
      free(ar1_aln.cigar);
      std::string rx1;
      if (ubam.GetZTag("RX", rx1)) {
        bam1.AddZTag("RX", rx1);
      }
      if (ar1_aln.XA)
        bam1.AddZTag("XA", std::string(ar1_aln.XA));

      // add num sub opt
      //b.AddIntTag("SB", ar.a[i].sub_n);
      bam1.AddIntTag("XS", ar1_sec_as);
      free(ar1.a);
      return bam1;
    } else {
      return ubam;
    }
  }
}

Segments PairEndBWA(const SeqLib::BWAWrapper& bwa, const Segments& segs, const int MIN_READL) {
  const int MAX_FRAG_DIST = 10000;
  Segments out;
  mem_alnreg_v ar1, ar2;
  mem_aln_t ar1_aln, ar2_aln;
  int ar1_sec_as = 0, ar2_sec_as = 0;

  std::string r1 = segs[0].Sequence();
  std::string r2 = segs[1].Sequence();
  uint8_t *qual1 = 0, *qual2 = 0;
  qual1 = bam_get_qual(segs[0].raw());
  qual2 = bam_get_qual(segs[1].raw());
  SeqLib::BamRecord bam1, bam2;
  //r1
  if (r1.length() < MIN_READL) {
    ar1.n = 0;
    memset(&ar1_aln, 0, sizeof(mem_aln_t));
    ar1_aln.rid = -1;
    ar1_aln.pos = -1;
    ar1_aln.flag |= 0x4;
  } else {
    ar1 = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                     r1.length(),r1.data());
    if (ar1.n > 0 ) {
      size_t f_pidx = 0;
      for (size_t idx = 0; idx < ar1.n; ++idx) {
        if (ar1.a[idx].secondary < 0) {
          f_pidx = idx;
        } else {
          ar1_sec_as = std::max(ar1_sec_as, ar1.a[idx].score);
        }
      }
      ar1_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, r1.length(), r1.c_str(), &ar1.a[f_pidx]);
    }
    free(ar1.a);
  }
  //r2
  if (r2.length() < MIN_READL) {
    ar2.n = 0;
    memset(&ar2_aln, 0, sizeof(mem_aln_t));
    ar2_aln.rid = -1;
    ar2_aln.pos = -1;
    ar2_aln.flag |= 0x4;
  } else {
    ar2 = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                     r2.length(),r2.data());

    if (ar2.n > 0 ) {
      size_t f_pidx = 0;
      for (size_t idx = 0; idx < ar2.n; ++idx) {
        if (ar2.a[idx].secondary < 0) {
          f_pidx = idx;
        } else {
          ar2_sec_as = std::max(ar2_sec_as, ar2.a[idx].score);
        }
      }
      ar2_aln = mem_reg2aln(bwa.GetMemOpt(), bwa.GetIndex()->bns, bwa.GetIndex()->pac, r2.length(), r2.c_str(), &ar2.a[f_pidx]);
    }
  }
  if (r1.length() >= MIN_READL) {
    free(ar1.a);
  }
  if (r2.length() >= MIN_READL) {
    free(ar2.a);
  }
  //set flags
  ar1_aln.flag |= 1;
  ar2_aln.flag |= 1;
  ar1_aln.flag |= 0x40;
  ar2_aln.flag |= 0x80;
  if (ar1.n == 0 && ar2.n==0) {
    out.push_back(segs[0]);
    out.push_back(segs[1]);
    return out;
  }  else {
    if (ar1.n == 0) {
      ar2_aln.flag &= 8;
      if (ar2_aln.is_rev)  {
        ar1_aln.flag |= 0x20;
      }
      bam1 = segs[0];
      bam2 = cpputil::BwaAlignment2BamRecord(ar2_aln, segs[1].Qname(), r2, qual2);
    }
    if (ar2.n == 0) {
      ar1_aln.flag &= 8;
      if (ar1_aln.is_rev)  {
        ar2_aln.flag |= 0x20;
      }
      bam1 = cpputil::BwaAlignment2BamRecord(ar1_aln, segs[0].Qname(), r1, qual1);
      bam2 = segs[1];
    }
  }
  if (ar1.n > 0 and ar2.n > 0) {
    if (ar1_aln.rid == ar2_aln.rid and (ar1_aln.is_rev ^ ar2_aln.is_rev)) {
      if (ar1_aln.is_rev and ar1_aln.pos >= ar2_aln.pos and ar1_aln.pos - ar2_aln.pos < MAX_FRAG_DIST) {
        ar1_aln.flag |= 2;
        ar2_aln.flag |= 2;
      } else if (ar2_aln.is_rev and ar2_aln.pos >= ar1_aln.pos and ar2_aln.pos - ar1_aln.pos < MAX_FRAG_DIST) {
        ar1_aln.flag |= 2;
        ar2_aln.flag |= 2;
      }
    }
    bam1 = cpputil::BwaAlignment2BamRecord(ar1_aln, segs[0].Qname(), r1, qual1);
    bam2 = cpputil::BwaAlignment2BamRecord(ar2_aln, segs[1].Qname(), r2, qual2);
  }
  //set other tag


  if (ar1.n > 0) {
    free(ar1_aln.cigar);
  }
  if (ar2.n > 0) {
    free(ar2_aln.cigar);
  }
  bam1.shared_pointer()->core.mtid = bam2.ChrID();
  bam1.shared_pointer()->core.mpos = bam2.Position();
  bam2.shared_pointer()->core.mtid = bam1.ChrID();
  bam2.shared_pointer()->core.mpos = bam1.Position();
  if (bam1.ProperPair()) {
    if (bam1.Position() < bam2.Position()) {
      bam1.shared_pointer()->core.isize = bam2.PositionEnd() - bam1.Position() + 1;
      bam2.shared_pointer()->core.isize = - bam2.PositionEnd() + bam1.Position() - 1;
    } else {
      bam1.shared_pointer()->core.isize = -bam1.PositionEnd() + bam2.Position() - 1;
      bam2.shared_pointer()->core.isize = bam1.PositionEnd() - bam2.Position() + 1;
    }
  }
  bam1.AddIntTag("XS", ar1_sec_as);
  bam2.AddIntTag("XS", ar2_sec_as);

  std::string rx1, rx2;
  if (segs[0].GetZTag("RX", rx1)) {
    bam1.AddZTag("RX", rx1);
  }
  if (segs[1].GetZTag("RX", rx2)) {
    bam2.AddZTag("RX", rx2);
  }
  if (r1.length() < MIN_READL) {
    SeqLib::Cigar c(std::to_string(r1.length()) + "S");
    bam1.SetCigar(c);
  }
  if (r2.length() < MIN_READL) {
    SeqLib::Cigar c(std::to_string(r2.length()) + "S");
    bam2.SetCigar(c);
  }
  out.push_back(bam1);
  out.push_back(bam2);
  if (bam1.CigarSize() > 0 and bam1.GetCigar().NumQueryConsumed() != r1.length()) {
    std::cerr << "invalid cigar " << bam1 << std::endl;
    out[0] = segs[0];
  }
  if (bam2.CigarSize() > 0 and bam2.GetCigar().NumQueryConsumed() != r2.length()) {
    std::cerr << "invalid cigar " << bam2 << std::endl;
    out[1] = segs[1];
  }
  return out;
}

}