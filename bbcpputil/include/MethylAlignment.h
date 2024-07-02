//
// Created by Ruolin Liu on 4/4/24.
//

#ifndef CODECSUITE_METHYLALIGNMENT_H
#define CODECSUITE_METHYLALIGNMENT_H

#include <string>
#include "BamRecordExt.h"
#include "Alignment.h"
namespace cpputil {

int GetProtectedStrandOrientation(const cpputil::Segments& segs) {
  int protected_strand_orientation =  0; // 1 for forward, 2 for reverse, 0 for unknown
  if (segs.size() == 2) {
    uint8_t* xm1 = bam_aux_get(segs[0].raw(),"XM");
    uint8_t* xm2 = bam_aux_get(segs[1].raw(),"XM");
    if (xm1 and !xm2) {
      protected_strand_orientation = int(segs[1].ReverseFlag()) + 1;
    } else if (!xm1 and xm2) {
      protected_strand_orientation = int(segs[0].ReverseFlag()) + 1;
    }
  } else {
    // for single strand read, check XM tag. If XM tag does not exist, it is a protected strand
    uint8_t* xm1 = bam_aux_get(segs[0].raw(),"XM");
    if (!xm1) {
      protected_strand_orientation = int(segs[0].ReverseFlag()) + 1;
    }
  }
  return protected_strand_orientation;
}

}

#endif //CODECSUITE_METHYLALIGNMENT_H
