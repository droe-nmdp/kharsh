package org.nmdp.ngs.kharsh

/*
 * Gibbs sampling for read origin on a fixed haplotype H1 and 
 * reference haplotype S. i.e., P(rj = -1 | H, S). It returns the probability
 * that the read is from haplotype h, as opposed to its complement h-bar.
 * Equation 6 in the HARSH paper. S is apparently not explictly used
 * in the sampling of the reads.
 * 
 * R is a vector of N indicator variables in {-1,1} where N is the number
 * of reads. -1 means the read aligns to the first haplotype, 1 means it
 * aligns to the second haplotype. 
 *
 * H is a two row matrix. Each row is a vector representing a predicted
 * haplotype. Each column is one of M indicator variables in {-1,1} where
 * M is the number of variants. -1 means the first variant and 1 means the 
 * second variant.
 *
 * hMarker holds the marker for H1 and H2 (e.g., -1 and 1) used in R.
 *
 * hIndex is an index into H (and hMarker) that indicates the haplotype to which the read is currently assigned
 *
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   Rows are the N reads(j), columns are the N variants(i).
 *
 * rIndex is the row index to X (and array index in R) indicating the sampled
 * read. 0-based index.
 *
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * @author Dave Roe
 * @version $Id: PFRScript.groovy 25301 2015-08-23 16:12:49Z droe $
 *
 */

class PFRScript {
    static err = System.err
    static int debugging = 1
        
    static Double PFR(int[][] H, int hIndex, int[][] X, int rIndex,
                      Double epsilon) {
        if(debugging <= 1) {
            err.println "PFR(hIndex=${hIndex}, rIndex=${rIndex})"
        }
        int otherhIndex = hIndex * -1 + 1  // either 0 or 1
        int[] hAssign = H[hIndex]
        int[] hOther = H[otherhIndex]

        // get the correct row (read)
        int[] xj = X[rIndex] // get the row vector for read at index rIndex
        // get the indexes of the variants assigned to this haplotype
        err.println "xj=${xj}" //droe
        err.println "hAssign=${hAssign}" //droe
        err.println "hOther=${hOther}" //droe
        int[] xRIndexes = [xj, hAssign].transpose().collect { it[0] == it[1] }.findIndexValues { it == true }
        // get the indexes of the variants assigned to the other haplotype
        int[] xROtherIndexes = [xj, hOther].transpose().collect { it[0] == it[1] }.findIndexValues { it == true }
        if(debugging <= 2) {
            err.println "xRIndexes=${xRIndexes}"
            err.println "xROtherIndexes=${xROtherIndexes}"
        }

        // the extent to which the variants assigned via this read (in X)
        // to this haplotype, match this haplotype (as indicated by hIndex)
        Double thetaSumhAssign =
            ThetaSumReadsScript.ThetaSumReads(xRIndexes, xj, hAssign, epsilon)
        // the extent to which the variants assigned via this read (in X)
        // to the other haplotype, do not match this haplotype
        Double nuSumhAssign =
            NuSumReadsScript.NuSumReads(xROtherIndexes, xj, hAssign, epsilon)
        // the extent to which the variants assigned via this read (in X)
        // to this haplotype, match the other haplotype
        // left off: the 'others' should be positive
        Double thetaSumhOther =
            ThetaSumReadsScript.ThetaSumReads(xRIndexes, xj, hOther, epsilon)
        // the extent to which the variants assigned via this read (in X)
        // to the other haplotype, do not match the other haplotype
        Double nuSumhOther =
            NuSumReadsScript.NuSumReads(xROtherIndexes, xj, hOther, epsilon)

        if(debugging <= 2) {
            err.println "PFR: thetaSumhAssign=${thetaSumhAssign}, nuSumhAssign=${nuSumhAssign}, " +
                "thetaSumhOther=${thetaSumhOther}, nuSumhOther=${nuSumhOther}"
        }

        
        Double numerator = Math.exp(thetaSumhAssign + nuSumhAssign)
        Double denominator = numerator + Math.exp(thetaSumhOther + nuSumhOther)

        if(debugging <= 2) { 
            err.println "PFR: ${numerator}/${denominator}"
        }
  
        Double prob = numerator / denominator
        if(debugging <= 2) { 
            err.println "PFR: return ${prob}"
        }
        return prob
    } // PFR
} // PFRScript