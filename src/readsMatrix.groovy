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
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
freqThreshold = 0.1 // ignore alleles less than this frequency
err = System.err

Map loadVariantMap(file) { 
  IOScript io = new IOScript()
  // position -> base -> representation (null, 1, -1)
  io.loadVariantMap(file)
}

class ReadState {
  Table<String, Integer, Integer> matrix1
  Table<String, Integer, Integer> matrixNeg1
  Integer previousRefPos
  String previousBase
}

processLine = 
{ Map variantMap, boolean translateVariants, Closure translateBase, ReadState state, String line
->
  def (name, flag, chrom, pos, base, qual, inRefPos, ref, op) = line.split('\t')
  if ( name =~ /READ_NAME/ ) {
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
  
  if (refPos < 0 ) { // todo (examine the meaning of this)
    state.previousRefPos = refPos
    state.previousBase = base
    return
  }
  
  if ((debugging <= 2) && ((base == '.') || (base == '*'))) {
    err.println "deletion at ref pos ${refPos}"
  }
  
  Integer variant = translateBase(refPos, base)//droe
  if ( variant != null ) {
    def putVariant = variant
    if ( translateVariants == false ) {
      putVariant = base
    }
    
    if ( variant == -1 ) {
      if ( debugging <= 2 ) { 
        err.println "adding: ${name}, ${refPos}, ${putVariant}"
      }
      
      state.matrixNeg1.put( name, refPos, putVariant )
    } else if ( variant == 1 ) {
      if ( debugging <= 2 ) { 
        err.println "adding: ${name}, ${refPos}, ${putVariant}"
      }
      state.matrix1.put( name, refPos, putVariant )
    } else {
      if ( debugging <= 2 ) { 
        err.println "warning unknown base: ${name}, ${refPos}, ${base}"
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
Integer processInputFile( file, process ) {
  Reader reader = new BufferedReader( file )
  reader.eachLine { line, index -> process line; index +1 }
}                                                                    

void main () {
  (inFile, outFile, vInFile, vOutFile, translateVariants) = handleArgs(args)
  
  err.println "loading input file ..."
  
  variantMap = loadVariantMap( vInFile )
  // variantMap: Gene index -> (map of base to marker)
  err.println "M: " + variantMap.size()

  Table<String, Integer, Integer> matrix1 = HashBasedTable.create()
  Table<String, Integer, Integer> matrixNeg1 = HashBasedTable.create()

  def state = new ReadState()
  state.matrix1 = matrix1
  state.matrixNeg1 = matrixNeg1
  
  readCount =
    processInputFile( inFile, processLine.curry( variantMap, translateVariants, translateBase.curry( vInFile == null, variantMap ), state ))
  
  err.println "done with reads"
  if ( debugging <= 3 ) { 
    err.println variantMap.keySet().size() + " positions in variantMap before pruning"
    err.println "matrixNeg1 positions=" + matrixNeg1.columnKeySet().size()
    err.println "matrix1 positions=" + matrix1.columnKeySet().size()
    if ( matrix1.columnKeySet().size() > 0 ) {
      // (todo) fix this: wrong arg type for ColumnKeySet.getAt
      ;//err.println "first matrix1 position=" + matrix1.columnKeySet()[(Integer)0]
    }
  }
  
  intersectSNPs = matrixNeg1.columnKeySet().intersect(matrix1.columnKeySet()).sort()

  err.println "I: " + intersectSNPs
  err.println "${intersectSNPs.size()} intersect SNPs"
  rows = matrix1.rowKeySet() + matrixNeg1.rowKeySet()
  rows.sort()
  
  err.println "R: " + rows
  err.println "${rows.size()} rows"
  Table<String, Integer, Integer> matrixCombined = HashBasedTable.create()
  intersectSNPs.each { col ->
    Map alleleCountMap = [:]
    Integer totalCount = 0
    rows.each( { row -> 
      value = matrixNeg1.get(row, col)
      if (value == null ) {
        value = matrix1.get(row, col)
      }
      
      if ( value != null ) {
        if ( debugging  <= 1 ) { 
          err.println "value=${value}"
        }
        matrixCombined.put(row, col, value)
        count = alleleCountMap[value]
        if ( count == null ) {
          alleleCountMap[value] = 1
        } else {
          alleleCountMap[value] = count+1
        }
        totalCount++
          }
               }    
      )
    err.println "T: ${totalCount}"
    // delete read variants if frequency is below a threshold
    if(vInFile != null) { // only restrict on the reads, not the reference alleles
      if(debugging  <= 3) { 
        err.println "removing alleles below ${freqThreshold} ..."
      }
      alleleCountMap.each { key, value->
        if((value/totalCount) < freqThreshold) {
          if(debugging  <= 1) { 
            err.println "removing col=${col}, key=${key}"
            err.println "value=${value}"
            err.println "totalCount=${totalCount}"
          }
          rows.each { row ->
            value2 = matrixCombined.get(row, col)
            if(value2 == key) {
              if(debugging <= 2) { 
                err.println "removing from matrixCombined ${row} ${col}"
              }
              matrixCombined.remove(row, col)
              removeList = []
              variantMap[col].each { key3, value3 ->
                if(debugging <= 2) { 
                  err.println "removing variantMap ${key3} (${value3}) for ${col}"
                  err.println variantMap[col].keySet()
                }
                removeList.add(key3)
              }
              removeList.each { key3 ->
                variantMap[col].remove(key3)
              }
            } //if the value matches
          } // each read
          if(debugging  <= 3) { 
            err.println "done removing alleles below threshold"
          }
        } // if threshold test
      } // each key (-1, 1, etc)
      if(variantMap[col].size() < 2) { // if only one left
        rows.each { row ->
          matrixCombined.remove(row, col)
        }
        variantMap.remove(col)
      }
    } // if reads run
  } // each snp/column
  matrixNeg1.clear()
  matrix1.clear()

  err.println "done. ${readCount} reads/sites, ${matrixCombined.rowKeySet().size()} reads in matrixCombined"
  err.println "${matrixCombined.columnKeySet().size()} snps in matrixCombined"

  err.println "outputing matrixCombined..."

  columns = matrixCombined.columnKeySet().sort()
  err.println columns.join(",")//todo
  // headers
  if(translateVariants == true) { 
    columns.each { column ->
      outFile.print "\t${column}"
    }
    outFile.println ""
  }
  // content
  rows = matrixCombined.rowKeySet().sort()
  rows.each { row ->
    if(translateVariants == true) {
      outFile.print "${row}"
    }
    columns.each { column ->
      value = matrixCombined.get(row, column)
      if(value == null) {
        value = ""
      }
      if(translateVariants == true) {
        outFile.print "\t"
      }
      outFile.print "${value}"
    }
    outFile.println ""
  }
  outFile.close()
  err.println "done"

  if(vOutFile != null) { 
    outputVariantMap(vOutFile, variantMap, matrixCombined)
  }
  //vOutFile.close()
}

/*
 * Returns a -1 or 1 representation of the read base. For now, the first
 * observed observation gets a -1 and the second gets a 1.
 * 
 * @param changeVariants if non-null, don't change variants; if null, change variants
 * @param Map variantMap position -> variant -> indicator
 * @param Integer pos the position of the variant
 * @param base the observed variant
 * @return a Integer in this order per position: -1, 1, 2, 3, ... or null
 *         for a non-variant position
 */
translateBase = 
{ boolean changeVariants, Map variantMap, Integer pos, String base
->

  if ( debugging <= 1) {
    err.println "translateBase(pos=${pos}, base=${base})"
  }
  
  Integer ret = null // not mapped or third+ variant
  baseMap = variantMap[pos] // map: variant -> indicator
  if ( baseMap != null ) {
    removeKeys = []
    addKeyVals = [:]
    baseMap.each { baseI, indicatorI ->
      if ( debugging <= 2 ) {
        err.println "translateBase: base=${base}, baseI=${baseI}, indicatorI=${indicatorI}"
      }
      if (base.startsWith(baseI) ) {
        ret = indicatorI
        if(base.length() > baseI.length()) { 
          if(debugging <= 2) {
            err.println "${pos}: adding ${base}, removing ${baseI}"
          }
          addKeyVals[base] = ret
          removeKeys.add(baseI)
        }
      }
    } // each base at this position

    removeKeys.each { key ->
      baseMap.remove(key)
    }
    addKeyVals.each { key, val ->
      baseMap[key] = val
    }
    
    if (ret != null) {
      if(debugging  <= 1) {
        err.println "translateBase: return ${ret}"
      }
      return ret
    }
    
    maxValue = baseMap.values().max()
    if(maxValue == null) {
      ret = -1
    } else if(maxValue == -1) {
      ret = 1
    } else {
      ret  = maxValue + 1
    }
    if ( changeVariants ) {
      variantMap[pos][base] = ret
    }
  } else { // first variant at a position
    if ( changeVariants ) {// assign the first one to be -1
      ret = -1
      variantMap[pos] = [:]
      variantMap[pos][base] = ret
    } else { // this is not a relevant position
      ret = null
    }
  }
  
  if ( debugging  <= 1 ) {
    err.println "translateBase: return ${ret}"
  }
  
  ret
}

void outputVariantMap( FileWriter vout, Map variantMap, HashBasedTable matrix) {
  def sorted = matrix.columnKeySet().sort()
  def concat = { sep, pos ->
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

  sorted.inject( "", { sep, pos -> vout.print "${sep}${pos}"; "\t" } )
  vout.println ""
  sorted.inject( "", concat )
  
  vout.println ""
  vout.close()

} // outputVariantMap

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
  if (tmpVal != null) {
    varOutFileVal = new FileWriter(tmpVal)
  }
  
  translateVariantsVal = (Boolean) translateVariantsArg.getValue()
  if (translateVariantsVal == null) {
    translateVariantsVal = new Boolean(true)
  }
  
  [ new FileReader(inFileArg.getValue())
    , new FileWriter(outFileArg.getValue())
    , new FileReader(varInFileArg.getValue())
    , varOutFileVal
    , translateVariantsVal
  ]
} // handleArgs


main()


