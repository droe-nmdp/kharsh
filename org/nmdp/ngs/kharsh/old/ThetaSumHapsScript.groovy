package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of Thetas given one of more columns(SNPs) from
 * a predicted alignment matrix and a read.
 * Theta is an edge potential that favors the read r and haplotype h to be equal.
 * See equations 7 and 8 in the HARSH paper.
 *
 * R is a vector of N indicator variables in {-1,1} where N is the number
 * of reads. 1 means is aligns to the first first haplotype, -1 means it
 * aligns to the inverse haplotype. 
 *
 * X is a N by M matrix computed from the read alignment against the reference. 
 * X_ij =0 (not matched to the location), 1 (same allele as the
 *   reference), -1 (different allele as the reference). 
 *   Rows are the N reads(j), columns are the N SNPs(x).
 *
 * xi is a single column (SNP) from alignment matrix X
 * 
 * snp is the genotype value being summed over the haplotype
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * todo: remove the X option (?) used only by P()?
 * @author Dave Roe
 * @version $Id: ThetaSumHapsScript.groovy 22961 2015-06-06 01:12:40Z droe $
 *
 */

class ThetaSumHapsScript {
    static err = System.err
    static int debugging = 6
    
    static Double ThetaSumHaps(int[] R, int[] xi, int[][] X, int snp, int rjMultiplier,
                       Double epsilon) {

        if(debugging <= 1) { 
            err.println "ThetaSumHaps(snp=${snp}, rjMultiplier=${rjMultiplier})"
        }
        Double sum = 0.0;
        if((X == null) || (X.size() == 0)) { // calculate sum for one haplotype
            //err.println "xi=${xi}" //todo
            int[] snpIndexes = xi.findIndexValues { it == 1 }
            //err.println "snpIndexes=${snpIndexes}" //todo
            if(snpIndexes.size() > 0) {
                (0..snpIndexes.size()-1).each { i ->
                    int rj = R[snpIndexes[i]] * rjMultiplier
                    sum = ThetaScript.Theta(rj, snp, epsilon, sum)
                    if(debugging <= 1) {
                        err.println "ThetaSumHaps: sum=${sum}"
                    }
                }
            }
        } else { // calculate sum for both haplotypes (X)
            err.println "ERROR in ThetaSumHapsScript: implement this"
            System.exit(1)
        }

        if(debugging <= 1) {
            err.println String.sprintf("ThetaSumHaps: return %1.9f", sum)
        }

        return sum
    } // ThetaSumHaps
} // ThetaSumHapsScript
