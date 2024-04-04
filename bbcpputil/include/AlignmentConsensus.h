//
// Created by Ruolin Liu on 3/19/20.
//

#ifndef CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_
#define CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_

#include <string>
#include <vector>
#include <map>
#include <htslib/sam.h>
#include "BamRecordExt.h"
#include "Alignment.h"
#include "SeqLib/BWAWrapper.h"
#include "DNAUtils.h"

namespace cpputil {

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

inline bool is_bisulfite_converted(const std::string& seq, int MIN_READL) {
  int numC = cpputil::countChar(seq, 'C');
  int numG = cpputil::countChar(seq, 'G');
  if (seq.size() > MIN_READL and ((float) numC / seq.size() < 0.05 or (float) numG / seq.size() < 0.05 )) {
    return true;
  } else {
    return false;
  }
}

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
