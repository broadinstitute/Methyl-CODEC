# MCODECsuite
MCODECsuite is a software tool designed to process Methyl-CODEC data. Built on [CODECsuite](https://github.com/broadinstitute/CODECsuite), it includes functions for demultiplexing, adapter trimming, methyl-alignment, and single-fragment mutation calling (SFC), all implemented in C++14.

## Installation
Tested on Red Hat 7 and Ubuntu 18.04

prerequisite for C++ based programs. For snakemake workflow check out [here](./snakemake)
1. git
2. tested with gcc 5.2 and 7.3 with c++14 support
3. cmake 3.18.3 or above

First, recursive clone the repo and create a build directory which will holds the installion files and final executables.

`git clone --recursive git@github.com:broadinstitute/Methyl-CODEC.git && cd Methyl-CODEC && mkdir build`

Next, build the program with cmake.

`cd build && cmake .. && make`

After this, you should be able to see an executable named `codec` in the build folder you just created.

## Demultiplexing
MCODECsuite is expected to work with raw lane-level fastq.gz. This can be obtained from illumina [bcl2fastq](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html).
The first step is demultiplexing and it requires a sample sheet in csv format for each lane which looks the following.
Currently, we have used 3 barcodes. 

| SampleName | IndexBarcode1 | IndexBarcode2 |
|------------|---------------|---------------|
|Sample01|CTTGAACGGACTGTCCAC|CACCGAGCGTTAGACTAC|
|Sample02|GAGCCTACTCAGTCAACG|GTGTCGAACACTTGACGG|
|Sample03|AGCTTGTAAGGCAGGTTA|ACTGATCTTCAGCTGACT|


`codex demux -1 reads.r1.fastq.gz -2 reads.r2.fastq.gz -p sample_sheet.csv -o demux_outprefix `

Given the toy sample_sheet.csv and code this command will generate 
```
demux_outprefix.sample_A.1.fastq.gz, demux_outprefix.sample_A.2.fastq.gz
demux_outprefix.sample_B.1.fastq.gz, demux_outprefix.sample_B.2.fastq.gz
```

## Adapter trimming
After demultiplexing CODEC reads still contain in-situ sample barcode and adapter sequences. The next step is to trim 
these out since they could interfere alignment

`codec trim -1 demux_outprefix.sample_A.1.fastq.gz -2 demux_outprefix.sample_A.2.fastq.gz -o trim_outprefix -u 3 -U 3 -f 2 -t 2 -s sample_A` 

This tells the CODECsuite that first 3bp of a read is the UMI and to trim off the next two 2bp. 
The output files of the adapter trimming step looks like
```
trim_outprefix.sample_A.trim.bam
trim_outprefix.sample_A.trim.log
```
By default, single-end byproducts are also output to the `trim.bam`. To split the output use `-S/--split_bam_output`.

The bam file is standard uBam (unmapped bam) with additional tags
```
RX: UMI sequence from R1 and R2, concatenated by a hyphen
QX: UMI quality scores
bc: Index barcode sequence
s5: 5' adapter sequence (same as Index barcode)
q5: 5' adapter quality scores
s3: 3' adapter sequence (same as Index barcode of the mate)
q3: 3' adapter quality scores
sl: the rest of 3' adapter  sequence
ql: the rest of 3' adapter quality scores
```

## Methyl alignment
After adapter trimming. The Methyl-CODEC reads are aligned by `codec ms-align -b input.bam -o correct_product_prefix -r REF -q 30 -d 12 -R 'read_group_string' -i byproduct_prefix > output.log`
```
input.bam: the adapter trimmed uBAM.
REF: the reference genome
correct_product_prefix: output prefix for correct products
byproduct_prefix: output prefix for byproducts
```
The correct Methyl-CODEC reads are aligned by BWA-MEM and BWA-SW and have bismark-like tag: `XR`, `XG` and `XM` (see [Bismark](https://github.com/FelixKrueger/Bismark) for more inforamtion) for methylation extraction. The BAM works seaminglessly with `bismark_methylation_extractor`. There is an additional tag for Methyl-CODEC: `XC`, which has two possible values:
`XC=0` for product 1 and `XC=1` for product 2.

## Single fragment caller (SFC) and mutation rate computation

The SFC in MCODEC is similar to CODEC and shares the same interface, with minor modifications to account for methyl-conversion in one of the reads.
```
codec call -b input.mark_duplicated.bam -L highconfidentregions.bed  -r hg19.fa -n germline.bam -p lenient -o output

```

The output of the MCODEC SFC are
```
output.mutation_metrics.txt: includes SNV_rate, INDEL_rate and etc. 
output.variatns_called.txt: mutations from single fragments
output.context_count.txt: trinucleotide context and dinucleotide context counts
output.monomer_count.txt: mutation rate at monomer contexts with all four standard bases as reference (A,C,G,T) and all four bases + 5mC as alternative bases. 
```
MCODEDsuite output.variatns_called.txt has two additional columns compared to that in the CODECsuite: 

```
pstrand_orientation: 0 for product 1 (protected strand on the forward strand), 1 for product 2 (protected strand ont he reverse strand)
meth_char: Z,z, X,x, H,h, U,u (see [Bismark](https://github.com/FelixKrueger/Bismark) for more inforamtion)
```

The end-to-end pipeline is available at `snakemake/AdapV2/Snakefile`
