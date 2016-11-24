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
 *   Rows are the N reads(j), columns are the N SNPs(i).
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
 *
 * @author Dave Roe
 * @version $Id: PScript.groovy 22962 2015-06-08 00:49:45Z droe $
 *
 */

class PScript {
    static err = System.err
    static int debugging = 3
    
    static Double P(int[] R, int[][] H, int[] S1, int[] S2, ArrayList SI1,
                   ArrayList SI2, Integer numRefHaps, int[][] X, Double mu, 
                   Double epsilon, Double omega, Double rho, Double Z) {
        if(debugging <= 1) {
            err.println "P()"
        }

        int[] H1 = H[0]
        int[] H2 = H[1]
        Integer N = X.size()   // number of reads
        Double thetaSum = 0.0;
        // left off(todo): make like Alg 2
        (0..N-1).each { rIndex ->        // each of the N reads
            int[] xi = X[rIndex] // get the row vector for read at index rIndex
            thetaSum += PFRScript.PFR(H, X, R, rIndex, mu, epsilon)
        } // each read
        
        //start old
        int[] Rp = R.collect { it * 1 } 
        
        // H1 and S1
        int[] xi = []
        int rj = 0
        Double thetaSumH1 = ThetaSumReadsScript.ThetaSumReads(xi, X, rj, Rp, H1,
                                                             epsilon)
        Double nuSumH1 = NuSumReadsScript.NuSumReads(xi, X, rj, Rp, H1, epsilon)
        Double epsilonH1S1 = EpsilonSumScript.EpsilonSum(H1, S1, omega)
        Double tauSumS1 = TauSumScript.TauSum(SI1, numRefHaps, rho)

        // H2 and S2
        Double thetaSumH2 = ThetaSumReadsScript.ThetaSumReads(xi, X, rj, R, H2, 
                                                             epsilon)
        Double nuSumH2 = NuSumReadsScript.NuSumReads(xi, X, rj, R, H2, epsilon)
        Double epsilonH2S2 = EpsilonSumScript.EpsilonSum(H2, S2, omega)
        Double tauSumS2 = TauSumScript.TauSum(SI2, numRefHaps, rho)
    
        Double sum = thetaSumH1 + nuSumH1 + epsilonH1S1 + tauSumS1 + thetaSumH2 + 
                    nuSumH2 + epsilonH2S2 + tauSumS2
        if(debugging <= 2) {
            err.println "P: ${thetaSumH1}, ${nuSumH1}, ${epsilonH1S1}, ${tauSumS1}, " +
                "${thetaSumH2}, ${nuSumH2}, ${epsilonH2S2}, ${tauSumS2}"
            err.println String.sprintf("P: sum=%1.9f", sum)
        }
        Double prob = (1/Z) * (Math.exp(mu * sum))
                    
        if(debugging <= 1) {
            err.println String.sprintf("P: return %1.9f", prob)
        }
        return prob
    } // P
} // PScript
