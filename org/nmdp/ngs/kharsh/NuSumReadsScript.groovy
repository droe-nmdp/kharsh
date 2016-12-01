package org.nmdp.ngs.kharsh

/*
 * Return the probablilty that the read is assigned to the first 
 * predictioned haplotype (1), as opposed to the second haplotype (-1).
 *
 * Calculate the sum of Nus given a read and one of a pair of 
 * haplotype predictions.
 *
 * Nu is an edge potential that favors the read r and haplotype h to be opposite. 
 * (i.e., all opposite variants).
 * See equations 3, 5, and 6 in the HARSH paper and ThetaSumReadsScript.
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   columns are the M variants(i), rows are the N reads(j).
 * xj is one row (i.e., read) from X that spans one or more variants (i).
 * H is a predicted haplotype
 * S is a reference haplotype
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: NuSumReadsScript.groovy 25296 2015-08-22 16:48:22Z droe $
 */

class NuSumReadsScript {
    static err = System.err
    static int debugging = 6

    static Double NuSumReads(int[] indexes, int[] xj, int[] H,
                             Double epsilon) {
        if(debugging <= 1) {
            err.println "NuSumReads(xj=${xj})"
        }
        Double sum = 0.0

        indexes.each { j ->
            int xjVal = xj[j]
            int hVal = H[j]
            sum += NuScript.Nu(xjVal, hVal, epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("NuSumReads: return %1.9f", sum)
        }
        return sum        
    } // NuSumReads
} // NuSumReadsScript