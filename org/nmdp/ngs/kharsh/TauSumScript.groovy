package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of Taus given two sites on a known/reference
 * haplotype.
 *
 * Tau models the "transition probability in haplotype copying model (Li
 * and Stephens, 2003)".
 *
 * S is a vector of indexes into the reference haplotypes for each snp.
 *
 * N is # of reference haplotypes
 *
 * rho is the population recombination rate
 *
 * @author Dave Roe
 * @version $Id: TauSumScript.groovy 25437 2015-08-29 13:02:06Z droe $
 */

class TauSumScript {
    static err = System.err
    static int debugging = 6
        
    static Double TauSum(ArrayList S, int N, Double rho) {
        Double sum = 0.0
        Double val = 0.0

        // -1 for 0-based indexes
        // -1 to not include the last snp (comparing i and i+1)
        int endIndex = S.size() - 2
        if(endIndex >= 0) { // -1 has special meaning in a range
            (0..endIndex).each { i ->
                int s1 = S.get(i)
                int s2 = S.get(i+1)
                if(s1 == s2) {
                    val = Math.exp(-rho/N) + (1-Math.exp(-rho/N))/N
                } else { 
                    val = (1-Math.exp(-rho/N))/N
                }
            } // each variant
        } // if at least 3 snps
        sum = sum + val

        return sum
    } // TauSum
} // TauSumScript