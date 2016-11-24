package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of Nus given one of more columns(SNPs) from
 * a predicted alignment matrix and a read.
 * See equations 7 and 8 in the HARSH paper.
 *
 * Nu is an edge potential that favors the read r and haplotype h to 
 * be opposite.
 *
 * R is a vector of N indicator variables in {-1,1} where N is the number
 * of reads. 1 means is aligns to the first first haplotype, -1 means it
 * aligns to the inverse haplotype. 
 *
 * X is a N by M matrix computed from the read alignment against the reference. 
 * X_ij =0 (not matched to the location), 1 (same allele as the
 *   reference), -1 (different allele as the reference). 
 *   Rows are the N reads(j), columns are the N SNPs(x).
 * If X is null, calculate sum for one snp (xi); if non-null, 
 * calculate sum for both haplotypes (X).
 *
 * xi is a single column (SNP) from alignment matrix X
 * 
 * snp is the genotype value being summed over the haplotype
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: NuSumHapsScript.groovy 22961 2015-06-06 01:12:40Z droe $
 */

class NuSumHapsScript {
    static err = System.err
    static int debugging = 6

    static Double NuSumHaps(int[] R, int[] xi, int[][] X, int snp,
                           int rjMultiplier, Double epsilon) {
        if(debugging <= 1) {
            err.println "NuSumHaps(snp=${snp})"
            err.println "NuSumHaps(xi=${xi})"//todo:remove
        }

        Double sum = 0.0
        if((X == null) || (X.size() == 0)) { // calculate sum for one snp (xi)
            int[] snpIndexes = xi.findIndexValues { it == -1 }
            if(snpIndexes.size() > 0) { 
                (0..snpIndexes.size()-1).each { i ->
                    int rj = R[snpIndexes[i]] * rjMultiplier
                    sum = NuScript.Nu(rj, snp, epsilon, sum)
                } // each xi snp with -1 as a value
            }
        } else { // calculate sum for all snps in X
            err.println "ERROR in ThetaSumHapsScript: implement this"
            System.exit(1)            
        } // X or xi
        
        if(debugging <= 1) {
            err.println "NuSumHaps: return ${sum}"
        }
        return sum
    } // NuSumHaps
} // NuSumHapsScript
