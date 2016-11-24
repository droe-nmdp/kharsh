package org.nmdp.ngs.kharsh

/*
 * Gibbs sampling for read origin on a fixed haplotype H1 and 
 * reference haplotype S. i.e., P(rj = -1 | H, S). It returns the probability
 * that the read is from haplotype h, as opposed to its complement h-bar.
 * Equation 5 in the HARSH paper
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
 * X is a N by M matrix computed from the read alignments against the reference. 
 * X_ij =0 (no variant), -1 (first variant), 1 (second variant). 
 *   Rows are the N reads(j), columns are the N variants(i).
 *
 * rIndex is the row index to X and a index into R indicating the sampled
 * read. 0-based index.
 *
 * mu represents the 'heat' of the model (todo: remove)
 * epsilon represents the sequencing error rate; must be > 0.0
 *
 * rho(j) = 
 *   numerator: prob. that variant 1 is on haplotype -1 and variant -1 is not on haplotype -1
 *   denominator: numerator + prob. that variant 1 is on haplotype 1 and variant -1 is no on haplotype 1
 *
 * @author Dave Roe
 * @version $Id: PFRScript.groovy 25301 2015-08-23 16:12:49Z droe $
 *
 */

class PFRScript {
    static err = System.err
    static int debugging = 1
        
    static Double PFR(int[][] H, int[][] X, int rIndex, Double mu,
                      Double epsilon) {
        if(debugging <= 1) {
            err.println "PFR(rIndex=${rIndex})"
        }

        int[] H1 = H[0]
        int[] H2 = H[1]
        int[] xi = X[rIndex] // get the row vector for read at index rIndex
        int rj = 1
        // todo: look into this: I think this is wrong
        int[][] nX = [][] // can't be null
        int[] nJ = [] // can't be null
        /*
         * Numerator: loop over all reads at that position. If assigned to
         * the first haplotype check theta, if assigned tot he other haplotype,
         * check nu.
         * Denominator: same except add the other haplotype
         */

        Double thetaSumH1 = ThetaSumReadsScript.ThetaSumReads(xi, H1, epsilon)

        Double nuSumH1 = NuSumReadsScript.NuSumReads(xi, nX, rj, nJ, H1,
                                                    epsilon)
        Double thetaSumH2 = ThetaSumReadsScript.ThetaSumReads(xi, H2, epsilon)
        Double nuSumH2 = NuSumReadsScript.NuSumReads(xi, nX, rj, nJ, H2,
                                                    epsilon)
        //TODO: these are all 0.0?

        if(debugging <= 2) {
            err.println "PFR: thetaSumH1=${thetaSumH1}, nuSumH1=${nuSumH1}, " +
                "thetaSumH2=${thetaSumH2}, nuSumH2=${nuSumH2}"
        }
        Double numerator = Math.exp(thetaSumH1 + nuSumH1)
        Double denominator = numerator + Math.exp(thetaSumH2 + nuSumH2)

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