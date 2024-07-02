//
// Created by Ruolin Liu on 3/19/20.
//

#ifndef CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_
#define CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_

#include <string>
#include <vector>
#include <map>
#include <htslib/sam.h>
#include "seqan/align.h"
#include "seqan_init_score.h"
#include "BamRecordExt.h"
#include "Alignment.h"
#include "SeqLib/BWAWrapper.h"
#include "DNAUtils.h"

namespace cpputil {

typedef int TValue;
typedef Score<TValue, ScoreMatrix<Dna5, Default> > TScoringScheme;

inline void find_insert_(const SeqLib::Cigar &cigar, int left_cursor, std::map<int, int> &ins_len) {
  for (auto it = cigar.begin(); it != cigar.end(); ++it) {
    if (it->Type() == 'H') {
      continue;
    } else if (it->Type() == 'I') {
      ins_len[left_cursor] = std::max(ins_len[left_cursor], (int) it->Length());
    } else {
      left_cursor += it->Length();
    }
  }
}

inline bool is_bisulfite_converted(const std::string& seq, int MIN_READL, float MAX_G_RATE = 0.05) {
  int numG = cpputil::countChar(seq, 'G');
  if (seq.size() > MIN_READL and ((float)  numG / seq.size() < MAX_G_RATE )) {
    return true;
  } else {
    return false;
  }
}

inline int IsCorrectMSPairedReads1(const Segments& segs) {
  /*
   * return 0: for incorrect MS reads
   * 1: correct MS reads with first read is the protected strand
   * 2: correct MS reads with second read is the protected strand
   */
  if (segs.size() != 2) {
    throw std::runtime_error("IsCorrectMSPairedReads1: segs.size() != 2");
  }
  bool a = is_bisulfite_converted(segs[0].Sequence(), 15);
  bool b = is_bisulfite_converted(segs[1].Sequence(), 15);
  if (a ^ b) {
    if (not a) return 1;
    else return 2;
  } else {
    return 0;
  }
}

int IsCorrectMSPairedReads2(const Segments& segs, const TScoringScheme& ss, const int min_ol_len);

std::string GetConsensusTemplate(const Segments& segs, int32_t& ref_most_left);

std::pair<std::string, std::string>
    GetGappedSeqAndQual(const SeqLib::BamRecord &r, const int start, const std::string& consensus_template);

std::string MergePairSeq(const Segments &segs, const std::vector<std::string>& seqs, bool trim_overhang);
std::string MergePair(const Segments &segs, bool trim_overhang);
std::pair<std::vector<std::string>, std::vector<std::string>> GetPairPileup(const Segments &segs);

std::string CallingMetC(const SeqLib::RefGenome& ref,const SeqLib::BamHeader& bamheader, const Segments &segs, bool trim_overhang, int qcutoff, int eof);

SeqLib::BamRecord SingleEndBWA(const SeqLib::BWAWrapper& bwa, const SeqLib::BamRecord& ubam, const int MIN_READL = 15);
Segments PairEndBWA(const SeqLib::BWAWrapper& bwa, const Segments& segs, const int MIN_READL = 15);
//std::pair<std::string, std::string> PairSeqConsensus(const Segments &seg, bool trim_overhang, int qcutoff);

}

#endif //CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_
