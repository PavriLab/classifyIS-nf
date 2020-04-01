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
     nextflow run dmalzl/classifyIS-nf --sitesA A_IS.bed --bamA conditionA.bam --labelA WT --sitesB B_IS.bed  --bamB conditionB.bam --labelB KD

     Options:
      --sitesA        BED file containing initiation sites for condition A (usually WT)
      --bamA          BAM file containing NS-seq reads for condition A
      --labelA        label to use for condition A (Default: A)

      --sitesB        BED file containing initiation sites for condition B (usually treated)
      --bamB          BAM file containing NS-seq reads for condition B
      --labelB        label to use for condition B (Default: B)

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

if (!params.sitesA) {
  exit 1, "--sitesA was not specified"
} else {
  if (!file(params.sitesA).exists()) {
    exit 1, "--sitesA was specified but file does not exist"
  }
}

if (!params.bamA) {
  exit 1, "--bamA was not specified"
} else {
  if (!file(params.bamA).exists()) {
    exit 1, "--bamA was specified but file does not exist"
  }
}

if (!params.sitesB) {
  exit 1, "--sitesB was not specified"
} else {
  if (!file(params.sitesB).exists()) {
    exit 1, "--sitesB was specified but file does not exist"
  }
}

if (!params.bamB) {
  exit 1, "--bamB was not specified"
} else {
  if (!file(params.bamB).exists()) {
    exit 1, "--bamB was specified but file does not exist"
  }
}

if (!file(params.outputDir).exists()) {
  file(params.outputDir).mkdir()
}

log.info ""
log.info " parameters"
log.info " ======================"
log.info " sitesA                   : ${params.sitesA}"
log.info " bamA                     : ${params.bamA}"
log.info " labelA                   : ${params.labelA}"
log.info " sitesB                   : ${params.sitesB}"
log.info " bamB                     : ${params.bamB}"
log.info " labelB                   : ${params.labelB}"
log.info " FC                       : ${params.FC}"
log.info " t                        : ${params.t}"
log.info " filePrefix               : ${params.filePrefix}"
log.info " outputDir                : ${params.outputDir}"
log.info " ======================"
log.info ""

mergeChannel = Channel
                  .fromList([[file(params.sitesA), file(params.sitesB)])

quantificationChannel = Channel
                            .fromList([[file(params.bamA), params.labelA,
                                        file(params.bamB), params.labelB]])

process mergeSites {

  tag { params.filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*_IS.bed"

  input:
  set file(sitesA), file(sitesB) from mergeChannel

  output:
  file("${params.filePrefix}_IS.bed") into resultsMergeSites

  shell:
  '''
  cat !{sitesA} !{sitesB} | \
  sort -k1,1 -k2,2n | \
  bedtools merge -c 4 -o collapse > !{params.filePrefix}_IS.bed
  '''
}

process quantifyReads {

  tag { params.filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.count.tsv"

  input:
  set file(bamA), val(labelA), file(bamB), val(labelB) from quantificationChannel
  file(mergedSites) from resultsMergeSites

  output:
  set val(labelA), val(labelB), file(mergedSites), file("${params.filePrefix}.count.tsv") into resultsQuantifyReads

  shell:
  '''
  multiBamSummary BED-file -b !{bamA} !{bamB} \
                           --BED !{mergeSites} \
                           -l !{labelA} !{labelB} \
                           -o !{params.filePrefix}.counts.npz \
                           --outRawCounts !{params.filePrefix}.counts.tsv \
                           -p 16 \
                           --ignoreDuplicates
  '''
}

process processQuantification {

  tag { params.filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.master.tsv"

  input:
  set val(labelA), val(labelB), file(mergedSites), file(quantificationFile) from resultsQuantifyReads

  output:
  set val(labelA), val(labelB), file("${params.filePrefix}.master.tsv") into resultsProcessQuantification

  shell:
  '''
  processCounts.py -ct !{quatificationFile} \
                   -b !{mergedSites} \
                   -wt !{labelA} \
                   -kd !{labelB} \
                   -FC !{params.FC} \
                   -t !{params.t} \
                   -o !{params.filePrefix}.master.tsv
  '''
}

process plotting {

  tag { params.filePrefix }

  publishDir  path: "${params.outputDir}",
              mode: "copy",
              overwrite: "true",
              pattern: "*.pdf"

  input:
  set val(labelA), val(labelB), file(masterTable) from resultsProcessQuantification

  output:
  file("${params.filePrefix}.pdf") into resultsPlotting

  shell:
  '''
  '''
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
