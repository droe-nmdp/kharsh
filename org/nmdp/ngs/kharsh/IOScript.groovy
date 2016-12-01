package org.nmdp.ngs.kharsh

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
    static Integer debugging = 1
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
     * loadVariantMap
     * Loads a file containing the variants into a Map
     *
     * @param vInFile a file with variant locations and each variants values
     * @return a Map: position -> base -> representation (null, 1, -1) 
     */
    Map loadVariantMap(FileReader vInFile) {
        Map variantMap = [:]  // position -> base -> representation (null, 1, -1)
        if(debugging <= 3) {
            err.println "loading load predefined variantMap ..."
        }
        def vReader = new BufferedReader(vInFile)
        def positions = vReader.readLine().tokenize()
        def variantLists = vReader.readLine().tokenize()
        for(int i=0; i <= positions.size(); i++) {
            //err.println "before continue, i=${i}"//todo
            // e.g., C=-1,G=1,.GG=2
            def pos = positions[i]
            if(debugging  <= 1) { 
                err.println "i=${i}, pos=${pos}"
            }
            if(pos == null) {
                //err.println "before continue, i=${i}"//todo
                continue
            }
            pos = pos.toInteger()
            def variantList = variantLists[i]
            if(debugging  <= 1) { 
                err.println "variantList=${variantList}"
            }
            variantList.tokenize(',').each { variant ->
                if(debugging  <= 1) { 
                    err.println "variant=${variant}"
                }
                String var
                String ind
                (var, ind) = variant.split("=")
                Integer indInt = ind.toInteger()
                def posMap = variantMap[pos]
                if(posMap == null) {
                    variantMap[pos] = [(var):(indInt)]
                } else {
                    posMap[var] = indInt
                }
            }
        } // each position
        err.println variantMap.keySet().size() + " positions in initial variantMap"        
        return variantMap
    } // loadVariantMap

    /* 
     * Converts a matrix generated from readsMatrix.groovy into the 
     * .map file needed for HARSH.
     *
     * Input is two matrices (one for the person, one for the reference alleles)
     * Output is written to outputFileName.
     */
    void makeHarshVcf(String fileName, String alleleFileName,
                      FileReader vInFile, String outputFileName) {
        if(debugging <= 1) {
            err.println "makeHarshVcf(fileName=${fileName}, alleleFileName=${alleleFileName}, vInFile=${vInFile}, outputFileName=${outputFileName})"
        }
        // position -> base -> representation (null, 1, -1)
        Map variantMap = loadVariantMap(vInFile)
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

       PrintWriter outFile = new File(outputFileName).newPrintWriter()
       // output header
       outFile.println "##fileformat=VCFv4.2"
       outFile.println "#CHROM\tPOS\tID\tREF\tALT"
       //todo pass in the reference name from fasta
       String chromosome = "2DL4*0010201" 
       commonSnpIndexes.eachWithIndex { snpIndex, listIndex ->
           def ref
           def alt
           def vm = variantMap[snpIndex]
           (ref, alt) = vm.entrySet()
           outFile.println "${chromosome}\t${snpIndex}\tSNP\t${ref.getKey()}\t${alt.getKey()}"
       } // each common snp

       outFile.close()
       if(debugging <= 1) {
           err.println "makeHarshVcf: return"
       }
    } // makeHarshVcf

} // IOScript
