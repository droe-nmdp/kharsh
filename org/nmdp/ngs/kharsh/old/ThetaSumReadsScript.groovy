package org.nmdp.ngs.kharsh

/*
 * Return values based on whether or not the read should be assigned 
 * to one haplotype or another (-1/1)..
 *
 * Calculate the sum of Thetas given a read and one or more variants from
 * an alignment matrix and a predicted haplotype.
 *
 * Theta is an edge potential that favors the read r and haplotype h to be equal 
 * (i.e., the same variant).
 * See equations 3, 5, and 6 in the HARSH paper. 
 * Equation 6. First, I assume that the probability is misprinted in the paper,
 * and the true equation is P(rj = 1|H, S).
 *   Numerator
 *     For all the loci with a certain read that that have variant 1 (Xi=1),
 *     sum the thetas for haplotype h1. For all the loci within the read with 
 *     the other variant (Xi=-1), sum the nus for haplotype 1.
 *  Denominator
 *    Numerator plus the same on the other haplotype.
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   Rows are the N reads(j), columns are the M variants(i).
 * xj is one row (i.e., read) from X that spans one or more variants (i).
 *
 * H is a vector of M SNP variables in {-1,1}. 0 (unobserved), -1 (first allele), 
 * 1 (second allele).  
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * TODO: check implications of passing in null for X and H (e.g., in PFR())
 *
 * @author Dave Roe
 * @version $Id: ThetaSumReadsScript.groovy 25295 2015-08-22 16:47:25Z droe $
 *
 */

class ThetaSumReadsScript {
    static err = System.err
    static int debugging = 6
    
    static Double ThetaSumReads(int[] xj, int[] H, Double epsilon) {
        if(debugging <= 1) {
            err.println "ThetaSumReads()"
        }
        Double sum = 0.0
        // find all the non-zero indexes in either list, and call Theta for each
        // 
        def varIndexes = xj.findIndexValues { it != 0 }
        varIndexes.addAll(H.findIndexValues { it != 0 })

        varIndexes.each { i ->
            int hi = H[i]
            sum = ThetaScript.Theta(rj, hi, epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("ThetaSumReads: return %1.9f", sum)
        }
        return sum        
    } // ThetaSumReads
    
} // ThetaSumReadsScript
