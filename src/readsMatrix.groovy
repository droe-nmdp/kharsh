#!/usr/bin/env groovy

/*
 * Generates a hit matrix of reads to variants given a vertical textual
 * representation of an alignment file as represented by the output of 
 * sam2tsv.
 *
 * Output matrix
 *  row headers: each column is a reference position starting with 1
 *  column headers: each row is a read/reference name (e.g. 2DL4*00101)
 *  cell: a '/' separated value containing 
 *        1) 1-based integer reference alleles
 *        2) the DNA variant
 *        3) the cigar operator (I,D,M)
 *    (e.g., 1/A/M, 2/AA/I, 3/AA/D)
 *
 * e.g., Step 1: make the matrix and alleles for the reference alleles
 *   ./src/readsMatrix.groovy -i tutorial/ref/2DL4_2/2DL4-alleles_reads.txt -o output/test1_2DL4_ref-matrix.txt
 * e.g., Step 2: make the matrix for the individual
 *   ./src/readsMatrix.groovy -i tutorial/ref/2DL4_2/2DL4-alleles_reads_00104.txt -o output/test1_2DL4_ref-matrix_00104.txt
 * 
 * The TSV file uses 0-based read indexes (READ_POS) and
 * 1-based reference positions (REF_POS).
 *
 * @see https://github.com/lindenb/jvarkit/wiki/SAM2Tsv
 * @see https://samtools.github.io/hts-specs/SAMv1.pdf
 *
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *  dsh-commandline.jar: http://www.dishevelled.org/commandline/
 *
 * @todo create path to output files if they don't exist
 * @todo: gunzip the input file(s) if required
 *
 * @author Dave Roe
 * @author Nick Wormley
 */

import org.nmdp.ngs.kharsh.*
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable
import java.io.*;
import org.dishevelled.commandline.*
import org.dishevelled.commandline.argument.*

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
freqThreshold = 0.1 // ignore alleles less than this frequency

err = System.err
IOScript io = new IOScript()

String removeRefAlleleFromMatrix(
    Table<String, Integer, Integer> varPosMatrix,
    Table<Integer, Integer, String> refVarPosMatrix,
    Table<String, Integer, Integer> inVarPosMatrix,
    String variant, Integer position) {
  if(debugging <= 2) {
    err.println "removeRefAlleleFromMatrix(variant=${variant}, position=${position})"
  }

  Integer refVar = varPosMatrix.get(variant, position)
  if(debugging <= 2) {
    err.println "removeRefAlleleFromMatrix: ${refVar}"
  }

  refVarPosMatrix.remove(refVar, position)
  ret = varPosMatrix.remove(variant, position)
  return ret
} // removeRefAlleleFromMatrix

/*
 * Returns the variant reference given a variant and position. Adds new
 * variants and references to the tables if necessary. Uses variants
 * previously defined in inVarPosMatrix when present.
 *
 * @param varPosMatrix Table of variants references per variant and position
 * @param refVariantMap Table of variants per variant reference and position
 * @param inVarPosMatrix Table of variants with previously assigned references; use for consistency
 * @return the variant reference
 */
Integer getRefAlleleFromMatrix(
    Table<String, Integer, Integer> varPosMatrix,
    Table<Integer, Integer, String> refVarPosMatrix,
    Table<String, Integer, Integer> inVarPosMatrix,
    String variant, Integer position) {
  if(debugging <= 2) {
    err.println "getRefAlleleFromMatrix(variant=${variant}, position=${position})"
  }
  
  // previously defined?
  Integer refVar = varPosMatrix.get(variant, position)
  if(refVar == null) { // defined in input matrix?
    if((inVarPosMatrix != null) &&
       (refVar = inVarPosMatrix.get(variant, position))) {
        refVarPosMatrix.put(refVar, position, variant)
        varPosMatrix.put(variant, position, refVar)
    }
  }

  // if new variant, set the new reference to max value + 1
  // could cache max value if this gets too slow
  if(refVar == null) { 
    Map column = varPosMatrix.column(position)
    //err.println "getRefAlleleFromMatrix: column=${column}"
    if((column != null) && (column.size() > 0)) {
      refVar = column.values().sort()[-1] // get max value
      refVar++
    } else {
      refVar = 1  
    }

    if(debugging <= 2) {
      err.println "getRefAlleleFromMatrix: adding ${position} ${variant} ${refVar}"
    }
    refVarPosMatrix.put(refVar, position, variant)
    varPosMatrix.put(variant, position, refVar)
  } // if reference variant doesn't exist

  if(debugging <= 2) {
    err.println "getRefAlleleFromMatrix: return ${refVar}"
  }
  return refVar
} // getRefAlleleFromMatrix

/**
 * Processes a single line from the input.
 *
 * 
 * @param inVarPosMatrix Table containing the previous variant reference; use for consistency
 * @return the last reference position (for use with insertions relative to reference)
 */
List processLine( 
    Table<String, Integer, String> refAllelePosMatrix,
    Table<String, Integer, Integer> varPosMatrix,
    Table<Integer, Integer, String> refVarPosMatrix,
    Table<String, Integer, Integer> inVarPosMatrix,
    String line,
    Integer prevRefPos,
    String prevVar) { 
  def (readName, flag, chrom, pos, base, qual, inRefPos, ref, op) = line.split('\t')
  if ( readName =~ /READ_NAME/ ) { // skip the header row
    return [1, ""]
  }
  if(debugging <= 1) {
    err.println "processLine(readName=${readName}, pos=${pos}, base=${base}, ref=${ref}, inRefPos=${inRefPos}, ref=${ref}, op=${op})"
  }

  if(prevVar.contains("/")) { 
    (prevPos, prevVar, prevType) = prevVar.split('/')
  }

  Integer retRefPos = prevRefPos
  /* 
   * Only consider deletions(D), insertions(I), or (mis)matches(M)
   * see page 5 of the SAM specification file
   * https://samtools.github.io/hts-specs/SAMv1.pdf
   * old
   *   ignore skipped(N), soft(S) or hard(H) clipped, or padded(P) regions
   *   op ~= /[NSHP]
   */
  if (!(op =~ /[DIM]/)) {
    return [retRefPos, ""]
  }

  Integer refPos // mark the starting point
  String variant = new String() // the full resulting variant
  Character variantType = null // M,I,D  mismatch, insertion, deletion
  // deletion relative to the reference will have a 'D' for its CIGAR op
  if(op == 'D') {
    // remove the old variant
    prevRefAllele = removeRefAlleleFromMatrix(varPosMatrix, refVarPosMatrix,
                                              inVarPosMatrix, prevVar, prevRefPos)
    if(prevRefAllele == null) {
      refPos = inRefPos.toInteger() // mark the starting point
      variant = ref
    } else {
      refPos = prevRefPos // use the original starting point
      variant = prevVar + ref
    }
    retRefPos = refPos
    variantType = 'D'
  } else if(op == 'I') { // insertion relative to reference
    // the variant has to start at the previous reference position
    refPos = prevRefPos
    retRefPos = prevRefPos
    variant = prevVar + base
    variantType = 'I'
    // remove the old variant
    removeRefAlleleFromMatrix(varPosMatrix, refVarPosMatrix,
                              inVarPosMatrix, prevVar, prevRefPos)
  } else if(op == 'M') { // match or mismatch
    refPos = inRefPos.toInteger() // mark the starting point
    retRefPos = refPos
    if(base != ref) {
      variant = base
      variantType = 'M'
    }
  }
  if(variant == "") { // base = reference
    return [retRefPos, base]
  }

  // find or create in varPosMatrix and refVarPosMatrix
  Integer refAllele = getRefAlleleFromMatrix(varPosMatrix, refVarPosMatrix,
                                             inVarPosMatrix, variant, refPos)
  variant = "${refAllele}/${variant}/${variantType}" // (e.g., 1/A/M, 2/AA/I, 3/AA/D)
  if((debugging <= 2) && (variantType != null)) {
    err.println "processLine: adding ${readName} ${refPos} ${variant}"
  }

  refAllelePosMatrix.put(readName, refPos, variant)
  return [retRefPos, variant]
} // processLine

/**
 * Runs the given process closure on each line of the given file.
 *
 * @process Closure to run on each line of input file
 * @return the total number of lines processed.
 */
Integer processInputFile(file,
                         Table<String, Integer, String> refAllelePosMatrix,
                         Table<String, Integer, Integer> varPosMatrix,
                         Table<Integer, Integer, String> refVarPosMatrix,
                         Table<String, Integer, Integer> inVarPosMatrix) {
  Reader reader = new BufferedReader( file )
  Integer previousPos = 1 // 1-based reference
  String previousVar = ""
  lineCount = 0
  reader.eachLine { line ->
    (previousPos, previousVar) = processLine(refAllelePosMatrix,
                                             varPosMatrix, refVarPosMatrix,
                                             inVarPosMatrix, line,
                                             previousPos, previousVar)
    lineCount++
  }
  return lineCount
} // processInputFile

/**
 * Removes the variants at the given position.
 */
removeVariant = {
  Map refVariantMap, Integer pos
  ->
    removeList = []
    refVariantMap[pos].each { key, value ->
      if ( debugging <= 2 ) {
        err.println "removing refVariantMap ${key3} (${value3}) for ${col}"
        err.println refVariantMap[pos].keySet()
      }
      removeList.add(key)
    }

    removeList.each { key -> refVariantMap[pos].remove(key) }
}

/**
 * This function modifies the given variant matrix by removing
 * variants that occur at a frequency below the required threshold.
 * @todo remove
 */
void removeInfrequentVariants(
  Map alleleCountMap
  , List readNames
  , Table<String, Integer, Integer> matrix
  , Integer pos
  , Integer totalCount
  , Closure removeVariant
) {

  if (debugging  <= 3) {
    err.println "removing alleles below ${freqThreshold} ..."
  }

  alleleCountMap.each { marker, count ->
    if ( (count/totalCount) >= freqThreshold ) {
      return
    }

    if ( debugging <= 1 ) {
      err.println "removing col=${col}, key=${key}"
      err.println "value=${value}"
      err.println "totalCount=${totalCount}"
    }

    readNames.each { readName ->
      value = matrix.get(readName, pos)
      if ( value != marker ) {
        if ( debugging <= 2 ) {
          err.println "removing from matrixCombined ${row} ${col}"
        }

        matrix.remove(readName, pos)
        removeVariant pos

      }
    }

    if (debugging  <= 3) {
      err.println "done removing alleles below threshold"
    }

  }

}

/*
 * Returns a -1 or 1 representation of the read base. For now, the first
 * observed observation gets a -1 and the second gets a 1.
 *
 * @param changeVariants if false, don't change variants in refVariantMap; if true, change variants in refVariantMap
 * @param Map refVariantMap position -> variant -> indicator
 * @param Integer pos the position of the variant
 * @param base the observed variant
 * @return a Integer in this order per position: -1, 1, 2, 3, ... or null
 *         for a non-variant position
 * @todo remove
 */
translateBase =
{ boolean changeVariants,
  Map refVariantMap,
  Integer pos,
  String base
->

  if (debugging <= 1) {
    err.println "translateBase(pos=${pos}, base=${base})"
  }

  baseMap = refVariantMap[pos] // map: variant -> indicator

  def maybeAddVariant = {
    if ( changeVariants ) {
      def map = refVariantMap[pos]
      if ( map == null ) {
        map = refVariantMap[pos] = [:]
      }

      max = map.values().max()
      if ( max == null ) {
        ret = -1
      } else if ( max == -1 ) {
        ret = 1
      } else {
        ret  = max + 1
      }

      refVariantMap[pos][base] = ret
    }

  }

  // First variant at a position
  if ( baseMap == null ) {
    return maybeAddVariant()
  }

  // Next variant at position

  // removeKeys = []
  //addKeyVals = [:]

  baseMap.each { baseI, indicatorI ->
    if ( debugging <= 2 ) {
      err.println "translateBase: base=${base}, baseI=${baseI}, indicatorI=${indicatorI}"
    }

    if ( base == baseI ) { //base.startsWith( baseI ) ) {
      ret = indicatorI
      if ( base.length() > baseI.length() ) {
        if ( debugging <= 2 ) {
          err.println "adding ${base}, removing ${baseI}"
        }

        //addKeyVals[base] = ret
        //removeKeys.add(baseI)
      }
    }
  } // each base at this position

  /*
  removeKeys.each { key ->
    baseMap.remove(key)
  }

  addKeyVals.each { key, val ->
    baseMap[key] = val
  }
  */
  if ( ret != null ) {
    if ( debugging <= 1) {
      err.println "translateBase: return ${ret}"
    }
    return ret
  }

  ret = maybeAddVariant()

  if ( debugging  <= 1 ) {
    err.println "translateBase: return ${ret}"
  }

  ret
}

/*
 * handleArgs
 *
 * @param args the command line arguments
 * @return list of command line arguments in a usable form; use -h for more info
 */
List handleArgs( String[] args ) {
  USAGE = "interp.groovy [args]"
  Switch help = new Switch("h", "help", "display help message")

  inFileArg = new FileArgument( "i"
                                , "input file from sam2tsv"
                                , "sam2tsv input"
                                , true
                              )
  outFileArg = new FileArgument( "o"
                                 , "output"
                                 ,"output matrix file"
                                 , true
                               )
  varInFileArg = new FileArgument( "p"
                                   , "input reference variants (from '-m' option)  "
                                   , "reference variants (input)"
                                   , false
                                 )

  ArgumentList arguments = new ArgumentList( help
                                             , inFileArg
                                             , outFileArg
                                             , varInFileArg
                                           )
  CommandLine commandLine = new CommandLine(args)
  try {
    CommandLineParser.parse(commandLine, arguments)
    if ( help.wasFound() ) {
      Usage.usage( USAGE, null, commandLine, arguments, System.out )
      System.exit(-2)
    }
  } catch (CommandLineParseException e) {
    Usage.usage( USAGE, e, commandLine, arguments, System.err )
    System.exit(-1)
  }

  [ new FileReader(inFileArg.getValue())
    , new FileWriter(outFileArg.getValue())
    , varInFileArg.getValue()
  ]
} // handleArgs

void main () {
  IOScript io = new IOScript()
  FileReader inFile
  File vInFile // reference variant file (optional; null if not present)
  FileWriter outFile
  (inFile, outFile, vInFile) = handleArgs(args)  
  err.println "loading input file ..."

  // hash variant references per read and position
  // Table<read, position, ref variant> (e.g., <2DL4*00101, 102, 1>)
  HashBasedTable<String, Integer, String> refAllelePosMatrix = HashBasedTable.create()
  // hash variants references per variant and position
  // Table<variant, position, ref variant> (e.g., <A, 102, 1>)
  Table<String, Integer, Integer> varPosMatrix = HashBasedTable.create()
  // hash variants per variant reference and position
  // Table<ref variant, position, variant> (e.g., <1, 102, A>)
  Table<Integer, Integer, String> refVarPosMatrix = HashBasedTable.create()

  // load a variant map; use for consistent variant references
  // input variant position matrix
  Table<String, Integer, Integer> inVarPosMatrix = null
  if(vInFile != null) {
    (refAllelePosMatrix, varPosMatrix, refVarPosMatrix) =
      io.loadVariantMatrix(new FileReader(vInFile))
  }

  readCount =
    processInputFile(
      inFile
      , refAllelePosMatrix
      , varPosMatrix
      , refVarPosMatrix
      , inVarPosMatrix
    )
  err.println "done with reads"

  if ( debugging <= 3 ) {
    err.println refAllelePosMatrix.rowKeySet().size() +
        " reads or reference alleles"
    err.println refAllelePosMatrix.columnKeySet().size() + " variant positions"
  }

  err.println "outputting ..."
  io.saveVariantMatrix(outFile, refAllelePosMatrix)

  err.println "done"
} // main

main()
