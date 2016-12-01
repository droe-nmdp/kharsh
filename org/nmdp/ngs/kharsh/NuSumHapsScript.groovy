package org.nmdp.ngs.kharsh

/*
 * Nu is an edge potential that favors the read r and haplotype h to 
 * be opposite. i.e., the extent to which the variants assigned to the 
 * other haplotype in X don't match this haplotype.
 * Calculate the sum of Nus given one of more columns(SNPs) from
 * a predicted alignment matrix and a read.
 * See equations 7 and 8 in the HARSH paper.
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   columns are the M variants(i), rows are the N reads(j).
 *
 * xi is one variant; multiple reads('xi's) cover the variant position; ~R
 * H is a predicted haplotype
 * S is a reference haplotype
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: NuSumHapsScript.groovy 22961 2015-06-06 01:12:40Z droe $
 */

class NuSumHapsScript {
    static err = System.err
    static int debugging = 6

    static Double NuSumHaps(int[] indexes, int[] xi, int[] H,
                               Double epsilon) {

        if(debugging <= 1) { 
            err.println "NuSumHaps(xi=${xi})"
        }
        Double sum = 0.0;

        indexes.each { i ->
            int xiVal = xi[i]
            int hVal = H[i]
            sum += NuScript.Nu(xiVal, hVal, epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("NuSumHaps: return %1.9f", sum)
        }

        return sum
    } // NuSumHaps
} // NuSumHapsScript
