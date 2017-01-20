#!/usr/bin/env groovy

/*
 * Generates a hit matrix of reads to SNPs given an alignment file. It will
 * output > 2 variants per position. The first two variants by occurance will get
 * -1 and 1.
 *
 * todo: update this order
 * e.g., Step 1: for the reads; ~5.5 minutes
 *   ./readsMatrix.groovy -i input/test1_2DL4.reads.txt -o output/test1_2DL4_matrix.txt -m output/test1_2DL4_variants.txt 2> output/test1_2DL4_matrix_errs.txt
 * e.g., Step 2: for the reference alleles; ~13 seconds
 *   ./readsMatrix.groovy -i input/test1-alleles_2DL4.reads.txt -o output/test1-alleles_2DL4_matrix.txt -p output/test1_2DL4_variants.txt -m output/test1-alleles_2DL4_variants.txt 2> output/test1-alleles_2DL4_matrix_errs.txt
 *   Currently, this step further restricts the short-read snps to overlap witht the reference alleles. This is usually good because it weeds out false snps.
 *
 * The TSV file uses 0-based read indexes (READ_POS) and
 * 1-based reference positions (REF_POS).
 *
 * Requires
 * CsvFileReader todo
 *  guava.jar: https://github.com/google/guava
 *  dsh-commandline.jar: http://www.dishevelled.org/commandline/
 *
 * export PATH=$HOME/repos/dev/bioinformatics/projects/ngs/scripts/groovy/:$PATH
 * @todo: if position is not in a loaded variantMap, then the two haps have to be homozygous; incorporate this in the code
 * @todo: on the reference allele run, rewrite the reads data after further cutting the snps
 * @todo: gunzip the input file(s) if required
 * @todo: move freqThreshold to be a command line arg
 * @todo: modularize
 * @todo: not working with indels
 *
 * @author Dave Roe
 * @version $Id: readsMatrix.groovy 25440 2015-08-29 13:09:07Z droe $
 */

import org.nmdp.ngs.kharsh.*
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable
import java.io.*;
import org.dishevelled.commandline.*
import org.dishevelled.commandline.argument.*
import org.nmdp.research.bio.util.CsvFileReader

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3
freqThreshold = 0.1 // ignore alleles less than this frequency
err = System.err

Map maybeLoadVariantMap(file) {
  if ( file == null ) {
    return [:]
  }

  reader = new FileReader(file)
  IOScript io = new IOScript()
  // position -> base -> representation (null, 1, -1)
  io.loadVariantMap(reader)
}

class ReadState {
  Table<String, Integer, Integer> matrix1
  Table<String, Integer, Integer> matrixNeg1
  Integer previousRefPos
  String previousBase
}

/**
 * This function processes a single line from the input.
 */
processLine =
  { Map variantMap,
    boolean translateVariants,
    Closure translateBase,
    ReadState state,
    String line
->
  def (readName, flag, chrom, pos, base, qual, inRefPos, ref, op) = line.split('\t')
  if ( readName =~ /READ_NAME/ ) {
    return
  }

  // So far, we are only looking at sequence matches, not insertions
  // or deletions.

  if ( op != "M" ) {
    return
  }

  // deletion relative to the reference will have a '.' for its value
  Integer refPos = null
  if ( inRefPos == '.' ) { // this rows is an insertion relative to the reference
    base = state.previousBase + base
    refPos = state.previousRefPos
    if ( debugging <= 2 ) {
      err.println "insertion ${base} at ref pos ${refPos}"
    }
  } else {
    refPos = inRefPos.toInteger()
  }

  if ( refPos < 0 ) { // todo (examine the meaning of this)
    state.previousRefPos = refPos
    state.previousBase = base
    return
  }

  if ((debugging <= 2) && ((base == '.') || (base == '*'))) {
    err.println "deletion at ref pos ${refPos}"
  }

  Integer variant = translateBase(refPos, base)
  if ( variant != null ) {
    def putVariant = variant

    if ( !translateVariants ) {
      putVariant = base
    }

    def mat = null
    if ( variant == -1 ) {
      mat = state.matrixNeg1
    } else if ( variant == 1 ) {
      mat = state.matrix1
    }

    if ( mat != null ) {
      if ( debugging <= 2 ) {
        err.println "adding: ${readName}, ${refPos}, ${putVariant}"
      }
      mat.put( readName, refPos, putVariant )
    } else {
      if ( debugging <= 2 ) {
        err.println "warning unknown base: ${readName}, ${refPos}, ${base}"
      }
    }
  }

  previousRefPos = refPos
  previousBase = base
}

/**
 * Runs the given process closure on each line of the given file.
 * Returns the total number of lines processed.
 */
Integer processInputFile( file, Closure process ) {
  Reader reader = new BufferedReader( file )
  reader.eachLine { line, index -> process line; index +1 }
}

/**
 * This function removes the variants at the given position.
 */
removeVariant = {
  Map variantMap, Integer pos
  ->
    removeList = []
    variantMap[pos].each { key, value ->
      if ( debugging <= 2 ) {
        err.println "removing variantMap ${key3} (${value3}) for ${col}"
        err.println variantMap[pos].keySet()
      }
      removeList.add(key)
    }

    removeList.each { key -> variantMap[pos].remove(key) }
}

/**
 * This function modifies the given variant matrix by removing
 * variants that occur at a frequency below the required threshold.
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

/**
 * This function removes variants from the given matrix
 * at the given position if the given variant map only contains
 * a single variant at that position.
 * The position is also removed from the variant map.
 */
void removeSingleVariants(
  Map variantMap
  , List readNames
  , Table<String, Integer, Integer> matrix
  , Integer pos
) {

  if ( variantMap[pos]?.size() == 1 ) { // if only one left
    readNames.each { readName -> mat.remove(readName, pos) }
    variantMap.remove(pos)
  }
}


/*
 * Returns a -1 or 1 representation of the read base. For now, the first
 * observed observation gets a -1 and the second gets a 1.
 *
 * @param changeVariants if false, don't change variants in variantMap; if true, change variants in variantMap
 * @param Map variantMap position -> variant -> indicator
 * @param Integer pos the position of the variant
 * @param base the observed variant
 * @return a Integer in this order per position: -1, 1, 2, 3, ... or null
 *         for a non-variant position
 */
translateBase =
{ boolean changeVariants,
  Map variantMap,
  Integer pos,
  String base
->

  if ( debugging <= 1) {
    err.println "translateBase(pos=${pos}, base=${base})"
  }

  baseMap = variantMap[pos] // map: variant -> indicator

  def maybeAddVariant = {
    if ( changeVariants ) {
      def map = variantMap[pos]
      if ( map == null ) {
        map = variantMap[pos] = [:]
      }

      max = map.values().max()
      if ( max == null ) {
        ret = -1
      } else if ( max == -1 ) {
        ret = 1
      } else {
        ret  = max + 1
      }

      variantMap[pos][base] = ret
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

void outputMatrix(
  FileWriter outFile
  , Table<String, Integer, Integer> mat
  , boolean translateVariants
) {

  positions = mat.columnKeySet().sort()

  // headers
  if ( translateVariants ) {
    positions.each { pos -> outFile.print "\t${pos}" }
    outFile.println ""
  }

  // content
  readNames = mat.rowKeySet().sort()
  readNames.each { readName ->
    if ( translateVariants ) {
      outFile.print "${readName}"
    }

    positions.each { pos ->
      value = mat.get(readName, pos)
      if ( value == null ) {
        value = ""
      }
      if ( translateVariants ) {
        outFile.print "\t"
      }
      outFile.print "${value}"
    }
    outFile.println ""
  }

  outFile.flush()
  err.println "done"
}

void outputVariantMap(
  FileWriter vout
  , Map variantMap
  , HashBasedTable matrix
) {

  def header = { sep, pos ->
    vout.print "${sep}${pos}"
    "\t"
  }

  def content = { sep, pos ->
    valueMap = variantMap[pos]
    vout.print "${sep}"
    vout.print valueMap.inject( "",
                                { ret, k, v ->
                                  if ( ret != "" ) {
                                    ret += ","
                                  }
                                  ret += "${k}=${v}"
                                }
                              )
    "\t"
  }

  def sorted = matrix.columnKeySet().sort()

  sorted.inject( "", header )
  vout.println ""

  sorted.inject( "", content )
  vout.println ""

  vout.flush()
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
                                , "input"
                                , "input file"
                                , true
                              )
  outFileArg = new FileArgument( "o"
                                 , "output"
                                 ,"output matrix file"
                                 , true
                               )
  varInFileArg = new FileArgument( "p"
                                   , "variantsIn"
                                   , "variants input file"
                                   , false
                                 )
  varOutFileArg = new FileArgument( "m"
                                    , "variantsOut"
                                    , "variant output file"
                                    , false
                                  )
  translateVariantsArg = new BooleanArgument( "t"
                                              , "translate"
                                              , "translate variants to -1, 1"
                                              , false
                                            )

  ArgumentList arguments = new ArgumentList( help
                                             , inFileArg
                                             , outFileArg
                                             , varInFileArg
                                             , varOutFileArg
                                             , translateVariantsArg
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

  FileWriter varOutFileVal
  tmpVal = varOutFileArg.getValue()
  if ( tmpVal != null ) {
    varOutFileVal = new FileWriter(tmpVal)
  }

  translateVariantsVal = (Boolean) translateVariantsArg.getValue()
  if (translateVariantsVal == null) {
    translateVariantsVal = new Boolean(true)
  }

  [ new FileReader(inFileArg.getValue())
    , new FileWriter(outFileArg.getValue())
    , varInFileArg.getValue()
    , varOutFileVal
    , translateVariantsVal
  ]
} // handleArgs


void main () {
  (inFile, outFile, vInFile, vOutFile, translateVariants) = handleArgs(args)

  err.println "loading input file ..."

  variantMap = maybeLoadVariantMap( vInFile )
  err.println "M: " + variantMap

  Table<String, Integer, Integer> matrix1 = HashBasedTable.create()
  Table<String, Integer, Integer> matrixNeg1 = HashBasedTable.create()

  def state = new ReadState()
  state.matrix1 = matrix1
  state.matrixNeg1 = matrixNeg1

  def changeVariants = vInFile == null

  readCount =
    processInputFile(
      inFile
      , processLine.curry(
        variantMap
        , translateVariants
        , translateBase.curry( changeVariants, variantMap )
        , state
      )
    )

  err.println "done with reads"
  if ( debugging <= 3 ) {
    err.println variantMap.keySet().size() + " positions in variantMap before pruning"
    err.println "matrixNeg1 positions=" + matrixNeg1.columnKeySet().size()
    err.println "matrix1 positions=" + matrix1.columnKeySet().size()
    if ( matrix1.columnKeySet().size() > 0 ) {
      err.println "first matrix1 position=" + matrix1.columnKeySet().first()
    }
  }

  intersectSNPs =
    matrixNeg1.columnKeySet().intersect(matrix1.columnKeySet()).sort()

  err.println "I: " + intersectSNPs
  err.println "${intersectSNPs.size()} intersect SNPs"

  readNames = matrix1.rowKeySet() + matrixNeg1.rowKeySet()
  readNames = readNames.sort()

  err.println "R: " + readNames
  err.println "${readNames.size()} reads"
  Table<String, Integer, Integer> matrixCombined = HashBasedTable.create()

  def incrementCountFor = { map, key ->
    count = map[key]
    if ( count == null ) {
      map[key] = 1
    } else {
      map[key] = count+1
    }
  }

  //  def removeVariant = removeVariant.curry( variantMap )

  intersectSNPs.each { pos ->
    Map alleleCountMap = [:]
    Integer totalVariantCount = 0

    readNames.each { readName ->
      value = matrixNeg1.get(readName, pos)
      if ( value == null ) {
        value = matrix1.get(readName, pos)
      }

      if ( value == null ) {
        return
      }

      if ( debugging  <= 1 ) {
        err.println "value=${value}"
      }

      matrixCombined.put( readName, pos, value )
      incrementCountFor( alleleCountMap, value )
      totalVariantCount++
    }

    err.println "T: ${totalVariantCount}"

    // delete read variants if frequency is below a threshold
    if ( vInFile != null) { // only restrict on the reads, not the reference alleles
      removeInfrequentVariants( alleleCountMap
                                , readNames
                                , matrixCombined
                                , pos
                                , totalVariantCount
                                , removeVariant
                              )
    }

    removeSingleVariants(variantMap, readNames, matrixCombined, pos)

  }

  err.println "done. ${readCount} reads/sites, ${matrixCombined.rowKeySet().size()} reads in matrixCombined"
  err.println "${matrixCombined.columnKeySet().size()} snps in matrixCombined"

  err.println "outputing matrixCombined..."

  outputMatrix( outFile, matrixCombined, translateVariants )

  if ( vOutFile != null ) {
    outputVariantMap( vOutFile, variantMap, matrixCombined )
  }

}

main()
