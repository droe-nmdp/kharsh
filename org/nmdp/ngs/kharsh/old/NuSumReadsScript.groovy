package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of Nus given one of more rows from an alignment
 * matrix and a predicted haplotype.
 * Nu is an edge potential that favors the read r and haplotype h to 
 * be opposite. See equation 3 of HARSH paper.
 * See equation 6 of HARSH paper.
 *
 * X is a N by M matrix computed from the read alignment against the reference. 
 * X_ij =0 (unobserved SNP), 1 (one allele), -1 (other allele). 
 *   Rows are the N reads(j), columns are the N SNPs(i).
 * xi is one row from X
 * Either X or xi needs to be 0. This indicates if all reads or one should
 * be summed.
 *
 * rj indicates if read is mapped to H1 (-1) or H2 (1), or not
 * at all (0). If X is 0, rj is a vector; otherwise rj is a number.
 *
 * H is a vector of M SNP variables in {-1,1} (one, other) where M
 * is the number of SNPs. 
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * TODO: check implications of passing in null for X and H (e.g., in PFR())
 *
 * @author Dave Roe
 * @version $Id: NuSumReadsScript.groovy 25296 2015-08-22 16:48:22Z droe $
 */

class NuSumReadsScript {
    static err = System.err
    static int debugging = 6

    static Double NuSumReads(int[] xj, int rj, int[] H, Double epsilon) {
        if(debugging <= 1) {
            err.println "NuSumReads()"
        }
        Double sum = 0.0
        // all the -1 variants
        def varIndexes = xj.findIndexValues { it == -1 }
        varIndexes.each { i ->
            int hi = H[i]
            sum = NuScript.Nu(1, hi, epsilon)
        }

        if(debugging <= 1) {
            err.println String.sprintf("NuSumReads: return %1.9f", sum)
        }
        return sum        
            
            
        
    } // NuSumReads
} // NuSumReadsScript