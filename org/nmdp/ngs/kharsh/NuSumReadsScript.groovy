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
 * @version $Id: NuSumReadsScript.groovy 25296 2015-08-22 16:48:22Z droe $
 */

class NuSumReadsScript {
    static err = System.err
    static int debugging = 6

    static Double NuSumReads(int[] xi, int[] H, int rj, Double epsilon) {
        if(debugging <= 1) {
            err.println "NuSumReads()"
        }
        Double sum = 0.0
        // find all the non-zero indexes
        def varIndexes = xi.findIndexValues { it != 0 }
        varIndexes.each { i ->
            int hi = H[i]
            int j = xi[i]
            sum = NuScript.Nu(rj, j, hi, Epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("NuSumReads: return %1.9f", sum)
        }
        return sum        
            
            
        
    } // NuSumReads
} // NuSumReadsScript