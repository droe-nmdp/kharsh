package org.nmdp.ngs.kharsh

/*
 * Sample haplotype H for fixed reads (R) and reference 
 * haplotype (S): P(H|R,S). See equations 7 and 8 in the HARSH paper.
 * i.e., the extent to which the variants assigned to the 
 * other haplotype in X don't match this haplotype.
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   columns are the M variants(i), rows are the N reads(j).
 *
 * indexes is the indexes of the reads assigned to this haplotype
 * xi is one variant; multiple reads('xi's) cover the variant position; ~R
 * H is a predicted haplotype
 * S is a reference haplotype
 * 
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: ThetaSumHapsScript.groovy 22961 2015-06-06 01:12:40Z droe $
 *
 */

class ThetaSumHapsScript {
    static err = System.err
    static int debugging = 6
    
    static Double ThetaSumHaps(int[] indexes, int[] xi, int[] H,
                               Double epsilon) {

        if(debugging <= 1) { 
            err.println "ThetaSumHaps(xi=${xi})"
        }
        Double sum = 0.0;

        indexes.each { i ->
            int xiVal = xi[i]
            int hVal = H[i]
            sum += ThetaScript.Theta(xiVal, hVal, epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("ThetaSumHaps: return %1.9f", sum)
        }

        return sum
    } // ThetaSumHaps
} // ThetaSumHapsScript
