package org.nmdp.ngs.kharsh

/*
 * Calulate the probability of the distribution between two haplotypes,
 * a read alignment, and a set of reference sequences. P(H,R,S;X)
 * See equation 3 of HARSH paper.
 *
 * R is a vector of N indicator variables in {-1,1} where N is the number
 * of reads. -1 means is aligns to the first first haplotype, -1 means it
 * aligns to the second haplotype. 
 *
 * H1 and H2 are vectors, each representing a predicted haplotype. 
 * Each column is one
 * of M indicator variables in {-1,1} where M is the number of SNPs. The
 * other haplotype is the inverse of H1.
 *
 * S1 and S2 are two vectors, each representing a reference haplotype.
 *
 * SI1 and SI2 are two vectors, each representing indexes of reference haplotypes.
 *
 * X is a N by M matrix computed from the read alignment against the reference. 
 * X_ij =0 (unobserved SNP), 1 (major allele), -1 (minor allele). 
 *   Rows are the reads(j), columns are the variants(i).
 *
 * mu represents the 'heat' of the model
 * epsilon represents the sequencing error rate; must be > 0.0
 * omega represents the haplotype copying 'error'
 * rho is the population recombination rate
 *
 * Z is a normalization constant to ensure the sum of P() for all R, H
 * equal 1.
 *
 * Return a probability, where mu and epsilon are two parameters. 
 * e.g, exp(200) = very large number
 *      exp(2) = 7
 *      exp(0.2) = 1.22
 *
 * @author Dave Roe
 * @version $Id: PScript.groovy 22962 2015-06-08 00:49:45Z droe $
 *
 */

class PScript {
    static err = System.err
    static int debugging = 1
    
    static Double P(int[][] H, int[] R, int[] S1, int[] S2, ArrayList SI1,
                   ArrayList SI2, Integer numRefHaps, int[][] X, Double mu, 
                   Double epsilon, Double omega, Double rho, Double Z) {
        if(debugging <= 1) {
            err.println "P(Z=${Z})"
        }

        int[] H1 = H[0]
        int[] H2 = H[1]
        int[] hMarkers = [-1, 1] // -1 = H1, 1 = H2
        Integer N = R.size()   // number of reads
        Double sum = 0.0

        // theta and nu for both haplotypes
        (0..N-1).each { rIndex ->        // each of the N reads
            // current haplotype assign for read
            int rHIndex = hMarkers.findIndexOf{ it == R[rIndex] }
            Double prob = PFRScript.PFR(H, rHIndex, X, rIndex, epsilon)
            sum += prob
        } // each read
        if(debugging <= 3) {
            err.println String.sprintf("P: after reads, sum=%1.9f", sum)
        }
        
        // epsilon for both haplotypes
        Double epsilonH1S1 = EpsilonSumScript.EpsilonSum(H1, S1, omega)
        Double epsilonH2S2 = EpsilonSumScript.EpsilonSum(H2, S2, omega)
        Double tauSumS1 = TauSumScript.TauSum(SI1, numRefHaps, rho)
        Double tauSumS2 = TauSumScript.TauSum(SI2, numRefHaps, rho)

        // the sum should be positive, substracted by epsilon sum,
        // and then added by tauSum
        // todo: start with matches from the start (left off)
        sum += epsilonH1S1 + tauSumS1 + epsilonH2S2 + tauSumS2
        Double muVal = Math.exp(mu * sum)
        if(debugging <= 3) {
            err.println "P: ${epsilonH1S1}, ${tauSumS1}, " +
                "${epsilonH2S2}, ${tauSumS2}"
            err.println String.sprintf("P: sum=%1.9f", sum)
            err.println String.sprintf("P: muVal=%1.9f", muVal)
        }

        Double prob = 1/Z * muVal
                    
        if(debugging <= 1) {
            err.println String.sprintf("P: return %1.9f", prob)
        }
        return prob
        System.exit(0) //todo(left off)
    } // P
} // PScript
