#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2020 Daniel Malzl
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

def helpMessage() {
    log.info"""
    ================================================================
     classifyIS-nf
    ================================================================
     DESCRIPTION
     Quantification of nasent strand sequencing reads in called initiation sites
     and classification of initiation sites according to this data

     Usage:
     nextflow run dmalzl/classifyIS-nf --sitesA IS_A.bed --sitesB IS_B.bed --bamA conditionA.bam --bamB conditionB.bam

     Options:
      --sitesA        BED file containing initiation sites for condition A (usually WT)
      --sitesB        BED file containing initiation sites for condition B (usually treated)
      --bamA          BAM file containing NS-seq reads for condition A
      --bamB          BAM file containing NS-seq reads for condition B

      --FC            log2(RPM) cutoff to use for assigning upregulation or downregulation for a given peak (Default: 0.585 [log2(1.5)])
      --t             log2(RPM) cutoff to use for specifying putative peak call threshold (Default: 2)

      --filePrefix    prefix to use for the output files
      --outputDir     directory to write results to (Default: results)

     Authors:
     Daniel Malzl (daniel.malzl@imp.ac.at)
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}
