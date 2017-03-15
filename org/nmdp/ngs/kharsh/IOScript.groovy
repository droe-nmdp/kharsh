package org.nmdp.ngs.kharsh

import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

/*
 * IO-related operations.
 *
 * @author Dave Roe
 * @version $Id: IOScript.groovy 25439 2015-08-29 13:04:07Z droe $
 * @todo switch the debugging levels
 */

import org.nmdp.research.bio.util.CsvFileReader
import java.io.*;

class IOScript {
    static Integer debugging = 3 // TRACE=1, DEBUG=2, INFO=3
    static err = System.err

    /*
     * Loads the S and X matrices. Constrains SNPs to the ones in both the
     * reference variants (S) and the individual (X); among other possible reasons,
     * this is needed due to the required reciprocity of the two haplotypes.
     *
     * @param fileName Person file. path&name of file containing the S matrix of known/reference
     *                 haplotypes. Each row is a reference sequence and each column
     *                 is a SNP value {-1,1}. Contains both row and column headers.
     * @param alleleFileName Reference alleles file. path&name of file containing the X matrix computed from
     *                       the read alignment against the reference. Contains both 
     *                       row and column headers.
     * @return an ArrayList containing:
     *         int[][] X
     *         String[] readNames
     *         int[][] S
     *         String[] alleleNames
     *         List<Integer> commonSnpIndexes
     * @todo modularize
     */
    ArrayList LoadMatrix(String fileName, String alleleFileName) {
        if(debugging <= 1) {
            err.println "LoadMatrix(fileName=${fileName}, alleleFileName=${alleleFileName})"
        }
        // get the common snps
        String header

        File personFile = new File(fileName)
        header = personFile.withReader { header = it.readLine() }
        List<Integer> snpIndexes = header.split('\t')[1..-1] // skip header column
        int numXReads = -1 // -1 for the header
        personFile.eachLine { numXReads++ }

        File alleleFile = new File(alleleFileName)
        header = alleleFile.withReader { header = it.readLine() }
        List<Integer> snpAlleleIndexes = header.split('\t')[1..-1] // skip header column 
        List<Integer> commonSnpIndexes = snpIndexes.intersect(snpAlleleIndexes).collect{ it as Integer}
        int numSReads = -1 // -1 for the header
        alleleFile.eachLine { numSReads++ }
        if(debugging <= 4) {
            err.println "numXReads=${numXReads}, numSReads=${numSReads}"
            err.println "${snpIndexes.size()} person snps and ${snpAlleleIndexes.size()} allele snps"
            err.println "LoadMatrix: ${commonSnpIndexes.size()} common variants"
            err.println "LoadMatrix: commonSnpIndexes=${commonSnpIndexes}"
        }
        
        // load X, the matrix of reads to reference
        err.println "loading X via ${fileName}..."
        ArrayList<String> readNames = new ArrayList<String>(numXReads)
        def person_csv = new CsvFileReader(fileName, true)
        int rowIndex = -1
        int[][] X = new int[numXReads][commonSnpIndexes.size()]
        person_csv.forEachRow { personRow ->
            rowIndex++;
            //err.println "personRow=${personRow}"//todo
            //err.println "read name? " + personRow[""] //todo
            // the first column (blank header value) contains the read name
            readNames.add(personRow[""])
            commonSnpIndexes.eachWithIndex { snpLocation, idx ->
                //err.println "snpLocation=${snpLocation}" //todo
                def cell = personRow[snpLocation.toString()]
                //err.println "cell=${cell}" //todo
                if((cell != null) && (cell != "")) {
                    if(debugging <= 2) { 
                        err.println "processMatrix: adding X[${rowIndex}][${idx}] = " +
                            cell.toInteger()
                    }
                    X[rowIndex][idx] = cell.toInteger()
                    
                } else {
                    X[rowIndex][idx] = 0
                }
            } // each of the common snp indexes
        } // each row in the person's matrix file
        //err.println "X[0][0]=${X[0][0]}"  // todo

        
        // load S, the matrix of known alleles to reference
        err.println "loading S via ${alleleFileName}..."
        ArrayList<String> alleleNames = new ArrayList<String>(numSReads)
        def ref_csv = new CsvFileReader(alleleFileName, true)
        rowIndex = -1
        int[][] S = new int[numSReads][commonSnpIndexes.size()]
        ref_csv.forEachRow { alleleRow ->
            rowIndex++;
            if(rowIndex >= numSReads) { // e.g., skip blank lines
                return
            }
            if(debugging <= 2) {
                err.println "alleleRow=${alleleRow}"
                err.println "read name? " + alleleRow[""]
            }
            alleleNames.add(alleleRow[""])
            commonSnpIndexes.eachWithIndex { snpLocation, idx ->
                //err.println "snpLocation=${snpLocation}" //todo(remove)
                def cell = alleleRow[snpLocation.toString()]
                //err.println "cell=${cell}" //todo(remove)
                if((cell != null) && (cell != "")) {
                    if(debugging <= 2) { 
                        err.println "processMatrix: adding S[${rowIndex}][${idx}] = " +
                            cell.toInteger()
                    }
                    S[rowIndex][idx] = cell.toInteger()
                    
                } else {
                    S[rowIndex][idx] = 0
                }
            } // each of the common snp indexes
            if(debugging <= 3)  {
                def alleleName = alleleRow[""]
                err.println "${alleleName}=${S[rowIndex]}"
            }
        } // each row in the reference allele matrix file
        err.println "done."
        //err.println "S[0][0]=${S[0][0]}"  // todo(remove)
        //err.println "S=${S}"  // todo(remove)

        if(debugging <= 1) {
            err.println "LoadMatrix: return"
        }

        return [X, readNames, S, alleleNames, commonSnpIndexes]
    } // LoadMatrix

    ArrayList processMatrix(List<Integer[]> rows, List<Integer> snpAlleleIndexes,
                            List<Integer> commonSnpIndexes) { 
        err.println "processMatrix: rows[0] size = ${rows[0].size()}" //todo
        err.println "${commonSnpIndexes.size()} commonSnpIndexes"//todo
        //err.println "commonSnpIndexes=" + commonSnpIndexes//todo
        ArrayList<String> readNames = new ArrayList<String>(rows.size())
        int[][] rows2 = new int[rows.size()][commonSnpIndexes.size()]
        // loop through each _person_ row
        int rows2RIdx = 0 // row index for rows2
        rows.eachWithIndex { row, rIdx ->
            if(rIdx == 0) { // first row is header
                return
            }
            //err.println "row size = ${row.size()}" //todo
            // loop through each _allele_ column and weed out ones not used
            int rows2CIdx = 0 // column index for rows2
            boolean usedThisRow = false
            /*
            err.println "processMatrix: row[1] = ${row[1]}"//todo
            err.println "commonSnpIndexes[0]=${commonSnpIndexes[0]}, snpAlleleIndexes[0]=${snpAlleleIndexes[0]}"//todo
            */
            row.eachWithIndex { cell, cIdx ->
                Integer referenceIndex = snpAlleleIndexes[cIdx-1].toInteger()
                if(debugging <= 2) {
                    err.println "processMatrix: snpAlleleIndexes[0]=${snpAlleleIndexes[0]}"
                    err.println "processMatrix: cell=${cell}, cIdx=${cIdx}, referenceIndex=${referenceIndex}"
                }
                Integer commonSnpIndex = commonSnpIndexes.findIndexOf { it == referenceIndex }
                /*if(referenceIndex == 9857) {  //todo: remove
                        System.exit(1)
                }
                if(cIdx == 12) { //todo: remove
                    System.exit(1)
                    }*/

                if((cIdx > 0) && ((commonSnpIndex == null) || (commonSnpIndex < 0))) {
                    if(debugging <= 2) { 
                        err.println "skipping snp ${referenceIndex}; not in common snps"
                    }
                    return
                }
                //err.println "cell=${cell}, cIdx=${cIdx}"//droe
                if(cIdx == 0) { // first column is the read name
                    /*
                    int underscoreIndex = cell.indexOf("_");
                    // clean up read names
                    if (underscoreIndex != -1) { 
                        cell = cell.substring(0, underscoreIndex);
                    }
                    err.println "adding ${cell} to readNames"//todo
                    */
                    readNames.add(cell)
                    return
                } else if((cell != null) && (cell != "")) { 
                    //todo(remove) rows2[rIdx-1][cIdx-1] = cell.toInteger()
//todo                    rows2[rows2RIdx][rows2CIdx++] = cell.toInteger()
                    if(debugging <= 2) { 
                        err.println "variant ${referenceIndex} is at index ${commonSnpIndex} in commonSnpIndexes"//todo
                        err.println "processMatrix: adding rows2[${rows2RIdx}][${commonSnpIndex}]=" +
                            cell.toInteger() + " for read ${readNames[-1]}"
                    }
                    rows2[rows2RIdx][commonSnpIndex] = cell.toInteger()
                    /*System.exit(1)
                    if(commonSnpIndex == 11) {  //todo: remove
                        System.exit(1)
                    }*/
                    usedThisRow = true
                } else {
                    //todo(remove) rows2[rIdx-1][cIdx-1] = 0
                    //err.println "rows2RIdx=${rows2RIdx}, rows2CIdx=${rows2CIdx}"//todo
                    rows2[rows2RIdx][commonSnpIndex] = 0
                    usedThisRow = true
                }
            } // each column in the allele matrix
            if(usedThisRow) {
                rows2RIdx++
            }
        } // each row
        //err.println "processMatrix: rows2[0]=${rows2[0]}" //droe
        return [rows2, readNames]
    } // processMatrix

    /*
     * loadVariantMatrix
     * Loads a file containing the variants into matrices.
     * The columns need to be ordered by position.
     *
     * @param vInFile a file with variant locations and each variants values
     * @return ArrayList of 3 Tables
     */
    ArrayList<Table> loadVariantMatrix(FileReader vInFile) {
        if(debugging <= 1) {
            err.println "loading variant matrix ..."
        }
        // hash variant references per read and position
        // Table<read, position, ref variant> (e.g., <2DL4*00101, 102, 1>)
        HashBasedTable<String, Integer, String> refAllelePosMatrix = HashBasedTable.create()
        // hash variants references per variant and position
        // Table<variant, position, ref variant> (e.g., <A, 102, 1>)
        Table<String, Integer, Integer> varPosMatrix = HashBasedTable.create()
        // hash variants per variant reference and position
        // Table<ref variant, position, variant> (e.g., <1, 102, A>)
        Table<Integer, Integer, String> refVarPosMatrix = HashBasedTable.create()
        def reader = new BufferedReader(vInFile)

        // first row are positions
        TreeSet<Integer> positions = new TreeSet()
        ArrayList header = reader.readLine().split('\t')
        header[1..-1].each { pos ->
            positions.add(pos.toInteger())
        }

        // non-header rows
        reader.eachLine { line ->
            if(!(line.contains('/'))) { // skip if no variants
                return
            }
            ArrayList cols = line.split('\t')
            String read = cols[0]
            if(debugging <= 2) { 
                err.println "loadVariantMatrix: line=${line}, read=${read}"
            }
            cols[1..-1].eachWithIndex { cell, i ->
                if((cell == null) || (cell == "")) {
                    return
                }
                String rv, var, type
                (rv, var, type) = cell.split('/')
                Integer refVar = rv.toInteger()
                Integer position = positions[i]
                if(debugging <= 2) { 
                    err.println "loadVariantMatrix: cell=${cell}"
                    err.println "loadVariantMatrix: rv=${rv}, var=${var}, type=${type}"
                }
                refVarPosMatrix.put(refVar, position, var)
                varPosMatrix.put(var, position, refVar)
                refAllelePosMatrix.put(read, position, rv)
            }
        } // each line

        return[refAllelePosMatrix, varPosMatrix, refVarPosMatrix]
    } // loadVariantMatrix

    /*
     * saveVariantMatrix
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
     * @param outFile FileWriter to direct output
     * @param refAllelePosMatrix Table<variant, position, ref variant> (e.g., <A, 102, 1>)
     */
    void saveVariantMatrix(FileWriter vout,
                           HashBasedTable<String, Integer, String> refAllelePosMatrix) {
        def header = { sep, pos ->
            vout.print "\t${pos}"
        }

        def rowContent = { seq, read ->
            vout.print "${read}" // row header
            Map columnMap = refAllelePosMatrix.row(read)
            if ( debugging <= 2 ) {
                err.println "saveVariantMatrix: columnMap=" + columnMap
            }
            refAllelePosMatrix.columnKeySet().sort().each { pos ->
                String val = columnMap[pos]
                val = val ? val : ""
                vout.print "\t${val}"
            }
            vout.println ""
        }
        
        def positions = refAllelePosMatrix.columnKeySet().sort()
        def reads = refAllelePosMatrix.rowKeySet().sort()

        // for the first column/header containing read/allele names
        positions.inject( "", header )
        vout.println ""

        reads.inject( "", rowContent )
        vout.println ""

        vout.flush()
        
    } //saveVariantMatrix

} // IOScript
