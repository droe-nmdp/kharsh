package org.nmdp.ngs.kharsh

/*
 * Return the probablilty that the read is assigned to the first 
 * predictioned haplotype (1), as opposed to the second haplotype (-1).
 *
 * Calculate the sum of Thetas given a read and one of a pair of 
 * haplotype predictions. Don't have to call it twice; since it is a binary
 * question, the probablilty is relative to the other haplotype.
 *
 * Theta is an edge potential that favors the read r and haplotype h to be equal 
 * (i.e., the same variants).
 * See equations 3, 5, and 6 in the HARSH paper. 
 * Equation 6. First, I assume that the probability is misprinted in the paper,
 * and the true equation is P(rj = 1|H, S).
 *   Numerator
 *     For all the loci with a certain read that that have variant 1 (Xi=1),
 *     sum the thetas for haplotype h1. For all the loci within the read with 
 *     the other variant (Xi=-1), sum the nus for haplotype h1.
 *  Denominator
 *    Numerator plus the same on the other haplotype (h2).
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   Rows are the N reads(i), columns are the M variants(j).
 * xj is one row (i.e., read) from X that spans one or more variants (j).
 *
 * h is one o two predicted haplotypes {-1,1}.
 *
 * rj is the current assignment of the read relative to haplotype {-1,1}.
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: ThetaSumReadsScript.groovy 25295 2015-08-22 16:47:25Z droe $
 *
 */

class ThetaSumReadsScript {
    static err = System.err
    static int debugging = 6
    
    static Double ThetaSumReads(int[] xi, int[] H, int rj, Double epsilon) {
        if(debugging <= 1) {
            err.println "ThetaSumReads()"
        }
        Double sum = 0.0
        // find all the non-zero indexes
        def varIndexes = xi.findIndexValues { it != 0 }

        varIndexes.each { i ->
            int hi = H[i]
            int j = xi[i]
            sum = ThetaScript.Theta(rj, j, hi, Epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("ThetaSumReads: return %1.9f", sum)
        }
        return sum        
    } // ThetaSumReads
    
} // ThetaSumReadsScript