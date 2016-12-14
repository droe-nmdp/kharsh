#!/usr/bin/env nextflow

/*
 * Full-gene interpretation for 2DL4 using KHARSH.
 *
 * Input are paired-end fastq files. The first file must contain 'read1'
 * and the second file must contain 'read2'. Generally takes < 5 minutes.
 *
 * First, 'source env.bash', then ./2dl4.nf. 0080101+00902
 *
 * Name like 008-009.read1.fasta and 008-009.read2.fasta.
 * 
 *
 * requires
 *   - bwa mem
 *   - https://github.com/lindenb/jvarkit/
 *   - readsMatrix.groovy
 *     e.g., export PATH=$PATH:$HOME/repos/dev/bioinformatics/projects/ngs/scripts/groovy/
 *   - export NXF_OPTS="-Xms1G -Xmx4G"
 *   - KHARSH
 *     e.g., export CLASSPATH=$HOME/src/stats/matlab/graph_model/kharsh_0.1.jar:$CLASSPATH
 *
 * @author Dave Roe
 * @version $Id: $
 * @todo: left off: have to do alleles alignment before second call to readsMatrix.groovy
 * @todo: add subscribes for notifications
 */
// Can see output via cat */*/*/*-results.txt.

// things that may change per run
Integer targetedRegionSize = (int)10000 // KIR and HLA class II; For HLA class I, set to 1000
Integer trimSize = targetedRegionSize / 1 // x% of reads

// kharsh parameters
Double mu =  1   // "heat" model parameter
Double epsilon = 0.01 // sequencing error rate
Double omega = 0.01   // haplotype copying 'error'
Double rho = 0.001    // the population recombination rate

//(droe: testing)Double convergeCutoff = 0.0001
Double convergeCutoff = 0.001
Integer arsStart = 3000
Integer arsEnd = 7500

// things less likely to change
params.experiment = "tutorial"
params.reference = "tutorial/ref/bowtie_indexes/2DL4_index"
ref = file("${params.reference}")
refFasta = file("tutorial/ref/2DL4.fasta")
refVCF = file("tutorial/ref/2DL4-alleles_reads.txt")
refVariants = file("tutorial/ref/2DL4-alleles_variants2.txt")//todo
refMatrix = file("tutorial/ref/2DL4-alleles_matrix2_3c.txt")//todo
readNameSuffix = ".read"
raw = "${params.experiment}/raw/2DL4/4chunks"
raw1 = "${raw}/**read1*{fastq,fq,fastq.gz,fq.gz}"
raw2 = "${raw}/**read2*{fastq,fq,fastq.gz,fq.gz}"
reads1 = Channel.fromPath(raw1).ifEmpty { error "cannot find any reads matching ${raw1}" }.map { path -> tuple(sample(path), path) }
reads2 = Channel.fromPath(raw2).ifEmpty { error "cannot find any reads matching ${raw2}" }.map { path -> tuple(sample(path), path) }

alignmentReadPairs = Channel.create()
readPairs = reads1.phase(reads2).map{ read1, read2 -> [ read1[0], read1[1], read2[1] ] }.tap(alignmentReadPairs)

/*
 * trimFastq
 *
 * Trim the FASTQ files for the sake of speed.
 * Just take the 'head's of each paired file: keep one read for every 10 bases
 * in the targeted region. This will give a depth of approximately 50 for both
 * KIR and HLA, assuming 250x2 paired-end reads.
 * 
 * @see targetedRegionSize
 * @param readPairs Two files representing paired FASTQ files
 * @return trimmedPairs Two shorter versions of the input files
 * @todo handle gzipped fastq files
*/
process gunzipFastq {
  input:
    set s, file(r1), file(r2) from readPairs
  output:
    set s, file{"${r1}.fastq"}, file{"${r2}.fastq"} into unzippedPairs

  """
  gunzip -f -c ${r1} > ${r1}.fastq
  gunzip -f -c ${r2} > ${r2}.fastq
  """
} // gunzipFastq

process trimFastq {
  input:
    set s, file(r1), file(r2) from unzippedPairs
  output:
    set s, file{"${r1}_trimmed.fastq"}, file{"${r2}_trimmed.fastq"} into trimmedPairs

  """
  cat ${r1} > ${r1}_trimmed.fastq
  cat ${r2} > ${r2}_trimmed.fastq
  #todo head -n ${trimSize/2} ${r1} > ${r1}_trimmed.fastq
  #todo head -n ${trimSize/2} ${r2} > ${r2}_trimmed.fastq
  """
} // trimFastq

/*
 * fastq2BAM
 *
 * Align two paired-end fastq files. Remove the unmapped reads, sort, and index.
 * e.g.,
 * input: 2DL4-03010101-150201_1.bwa.read1.fastq.gz, 2DL4-03010101-150201_1.bwa.read2.fastq.gz
 * output: 2DL4-03010101-150201_2DL4_mapped.bam, 2DL4-03010101-150201_2DL4_mapped.bam.bai
* todo: move the -p value to be system-dependent/configurable
*/
process fastq2BAM {
  input:
    set s, file(r1), file(r2) from trimmedPairs
  output:
    set s, file {"${s}_mapped.bam"} into personBAM

    """
  bwa mem -t2 ${refFasta} ${r1} ${r2} > ${s}.sam 2> ${s}_err.txt
  samtools view -b -S ${s}.sam > ${s}.bam
  samtools sort -o ${s}_sorted.bam ${s}.bam
  samtools view -b -F0x4 ${s}_sorted.bam > ${s}_mapped.bam
  samtools index ${s}_mapped.bam
  #nice gzip -f ${s}.sam
  #todo rm ${s}.sam
  """
//  bowtie2 --very-sensitive --end-to-end -x ${ref} -1 ${r1} -2 ${r2} > ${s}.sam 2> ${s}_err.txt
} // fastq2BAM

/*
 * bam2Variant
 *
 * Converts an alignment file to a VCF-like file.
 *
 */
process bam2Variant {
  input:
    set s, file(r) from personBAM
  output:
    set s, file {"${s}_reads.txt"} into personVCF

  """
  sam2tsv -r ${refFasta} ${r} > ${s}_reads.txt
  """
} // bam2Variant

/*
 * variants2PersonMatrix
 *
 * Generates a hit matrix of _person_ reads to variants given an alignment file.
 * Also output a list of the variants.
 *
 * @param personVCF
 * @return personMatrix
 * @return personVariant
 * @todo don't need personVariant output
 */
process variants2PersonMatrix {
  input:
    set s, file(r) from personVCF
  output:
    set s, file {"${s}_2DL4_matrix.txt"} into personMatrix
    //set s, file {"${s}_2DL4_variants.txt"} into personVariant

  """
  readsMatrix.groovy -i ${r} -o ${s}_2DL4_matrix.txt -m ${s}_2DL4_variants.txt -p ${refVariants}
  """
} // variants2PersonMatrix

/*
 * refVariants
 *
 * Generates a variant matrix containing the variants for the reference alleles.
 * They are specific to this individual.
 *
 * @param personMatrix a matrix of reference variants per read
 * @return refVariant
 * @todo remove this: not needed any more
process refVariants {
  input:
    set s, file(r) from personVariant
  output:
    set s, file {"${s}-alleles_2DL4_matrix.txt"} into refMatrix
    set s, file {"${s}-alleles_2DL4_variants.txt"} into refVariant

  """
  readsMatrix.groovy -i ${refVCF} -o ${s}-alleles_2DL4_matrix.txt -p ${r} -m ${s}-alleles_2DL4_variants.txt
  """
} // variants2RefMatrix
*/

/*
 * kharsh
 *
 * @param personMatrix matrix file of the person's reads and variants
 * @param refMatrix matrix file of the reference alleles and variants
 * @return outputs the results to results.txt and stdout
 *
 * e.g., kharsh.groovy -m 0.00001 -e 0.01 -o 0.001 -r 0.0001 -c 0.1 -f 800 -g 900 -p input/2DL4-03010101-150201_2DL4_matrix.txt -a input/2DL4-alleles_matrix.txt
 */
process kharsh {
  input:
    set s, file(pm) from personMatrix
  output:
    set s, file {"${s}-results.txt"} into results
    stdout result

  """
  kharsh.groovy -m ${mu} -e ${epsilon} -o ${omega} -r ${rho} -c ${convergeCutoff} -f ${arsStart} -g ${arsEnd} -p ${pm} -a ${refMatrix} > ${s}-results.txt
  cat ${s}-results.txt
  """
  
} // kharsh

result.subscribe {
    println ""
    println it
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  def ret = name.substring(start, name.indexOf(readNameSuffix))
  System.err.println "sample: return ${ret}"//droe
  return ret
} // sample
