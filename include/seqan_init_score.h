//
// Created by Ruolin Liu on 8/28/23.
//

#ifndef CODECSUITE_INCLUDE_SEQAN_INIT_SCORE_H_
#define CODECSUITE_INCLUDE_SEQAN_INIT_SCORE_H_
//![header]
/*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de
 ============================================================================
  Copyright (C) 2010

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.
 ============================================================================
  Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
 ============================================================================
  Demonstration on how to initialize a scoring matrix programatically with:

   - one of the built-in matrices, here BLOSUM30
   - arbitrary values
   - a new, built-in matrix.
 ==========================================================================*/
//![header]
//![includes]
#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>   // For printing strings.
#include <seqan/score.h>    // The module score.

using namespace seqan;
typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;
typedef seqan::Iterator<TRow>::Type TRowIterator;
//![includes]

//![user-defined-matrix]
// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for the DNA alphabet.  For this, we first create a new tag.
struct UserDefinedMatrix {};

// Then, we specialize the class ScoringMatrix_ for the Dna5 alphabet.
template <>
struct ScoringMatrixData_<int, Dna5, UserDefinedMatrix>
{
  enum
  {
    VALUE_SIZE = ValueSize<Dna5>::VALUE,
    TAB_SIZE = VALUE_SIZE * VALUE_SIZE
  };

  static int const other_mm = -4;
  static int const BS_mm = 0;
  static inline int const * getData()
  {
    // The user defined data table.  In this case, we use the data from BLOSUM-30.
    static int const _data[TAB_SIZE] =
        {
            1, other_mm, BS_mm, other_mm, other_mm,
            other_mm, 1, other_mm, BS_mm, other_mm,
            other_mm, other_mm, 1, other_mm, other_mm,
            other_mm, other_mm, other_mm, 1, other_mm,
            other_mm, other_mm, other_mm, other_mm, other_mm
        };
    return _data;
  }

};
}  // namespace seqan
//![user-defined-matrix]

//![show-scoring-matrix]
// Print a scoring scheme matrix to stdout.
template <typename TScoreValue, typename TSequenceValue, typename TSpec>
void showScoringMatrix(Score<TScoreValue, ScoreMatrix<TSequenceValue, TSpec> > const & scoringScheme)
{
  // Print top row.
  for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
    std::cout << "\t" << TSequenceValue(i);
  std::cout << std::endl;
  // Print each row.
  for (unsigned i = 0; i < ValueSize<TSequenceValue>::VALUE; ++i)
  {
    std::cout << TSequenceValue(i);
    for (unsigned j = 0; j < ValueSize<TSequenceValue>::VALUE; ++j)
    {
      std::cout << "\t" << score(scoringScheme, TSequenceValue(i), TSequenceValue(j));
    }
    std::cout << std::endl;
  }
}

std::tuple<int,int,int,int> print_AG_CT_mismatch(const TAlign& align) {
  int n_C_T_mm = 0;
  int n_A_G_mm = 0;
  int n_other_mm = 0;
  int n_match = 0;
  const auto& row1 = seqan::row(align, 0);
  const auto& row2 = seqan::row(align, 1);
  auto it1_beg = begin(row1);
  auto it1_end = end(row1);
  auto it2_beg = begin(row2);
  auto it2_end = end(row2);
  for (; it2_beg != it2_end; ++it1_beg,++it2_beg) {
    if (!isGap(it2_beg)) break;
  }
  for (; it1_end != it1_beg; --it2_end,--it1_end) {
    if (!isGap(it1_end)) break;
  }
  for (;it1_beg != it1_end && it2_beg != it2_end; ++it1_beg, ++it2_beg) {
    if (*it1_beg != *it2_beg) {
      if (*it1_beg == 'C' && *it2_beg == 'T') {
        n_C_T_mm++;
      } else if (*it1_beg == 'A' && *it2_beg == 'G') {
        n_A_G_mm++;
      } else {
        n_other_mm++;
      }
    } else {
      n_match++;
    }
  }
  return std::make_tuple(n_C_T_mm, n_A_G_mm, n_other_mm, n_match);
}


#endif //CODECSUITE_INCLUDE_SEQAN_INIT_SCORE_H_
