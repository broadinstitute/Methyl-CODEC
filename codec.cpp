//
// Created by Ruolin Liu on 5/10/21.
//

#include <iostream>
#include <string>
#include<cstring>
#include <getopt.h>
#include <cassert>
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0.0"
#endif

int codec_demux(int argc, char **argv);
int codec_trim(int argc, char **argv);
int codec_consensus(int argc, char **argv);
int codec_accuracy(int argc, char **argv);


int print_help()
{
  std::cout<< "---------------------------------------------------\n";
  std::cout<< "Program: codec (concatenating original duplex for error correction analysis suite)\n";
  std::cout<< "Version: " << PACKAGE_VERSION << std::endl;
  std::cout<< "Usage:   codec <command> [options]\n";
  std::cout<< "Common command:      demux                      demultiplex fastq file.\n";
  std::cout<< "                     trim                       trim CODEC adapter sequence.\n";
  std::cout<< "                     consensus                  merging overlapping paired ends.\n";
  std::cout<< "                     accuracy                   per-base alignment accuracy.\n";
  std::cout<< "---------------------------------------------------\n";
  std::cout<< "Optional command:    annotate_bam_with_umis     add UMI to bam.\n";
  std::cout<< "                     trimbam                    trim ends of fragment.\n";
  std::cout<< "---------------------------------------------------\n";
  std::cout<< "Contact: ruolin@broadinstitute.org. "
              "Copyright: bloodbiopsy@broadinstitute.org 2020-2021. \n";
  return 1;
}

int main(int argc, char *argv[]) {
  int ret;
  if (argc < 2) return print_help();
  else if (strcmp(argv[1], "demux") == 0) ret = codec_demux(argc-1, argv+1);
  else if (strcmp(argv[1], "trim") == 0) ret = codec_trim(argc-1, argv+1);
  else if (strcmp(argv[1], "consensus") == 0) ret = codec_consensus(argc-1, argv+1);
  else if (strcmp(argv[1], "accuracy") == 0) ret = codec_accuracy(argc-1, argv+1);
  else {
    std::cerr << "[codec] unrecongnized command " << argv[1] << std::endl;
    print_help();
    ret = 1;
  }
  return ret;
}