package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of epsilons given a predicted haplotype and a
 * reference haplotype.
 *
 * Epsilon is an edge potential representing 'haplotype copying', which
 * is motivated by the idea that the predicted haplotype is a mosaic of 
 * reference haplotypes with a small number of differences.
 * 
 * H is the predicted haplotype represented as a vector of M indicator
 * variables in {-1,1} where M is the number of SNPs. 
 *
 * S is the predicted haplotype represented as a vector of M indicator
 * variables in {-1,1} where M is the number of SNPs. 
 *
 * omega represents a novelty rate (i.e., the extent to which the predicted
 * haplotype is expected to differ from reference)
 *
 * @author Dave Roe
 * @version $Id: EpsilonSumScript.groovy 22961 2015-06-06 01:12:40Z droe $
 *
 */

class EpsilonSumScript {
    static err = System.err
    static int debugging = 6
        
    static Double EpsilonSum(int[] H, int[] S, Double omega) { 
        Double sum = 0.0
        if(debugging <= 1) {
            err.println "EpsilonSum()"
        }

        if(debugging <= 2) {
            err.println "EpsilonSum: H=${H}"
            err.println "EpsilonSum: S=${S}"
        }
        
        // compare the variants between H and S
        int numEqual = 0
        [H, S].transpose().collect { 
            if(it[0] == it[1]) {
                numEqual++
            }
        }
        int numNotEqual = H.length - numEqual

        if(debugging <= 2) {
            err.println "EpsilonSum: numEqual=${numEqual}, numNotEqual=${numNotEqual}"
        }

        sum = sum + Math.log(1-omega) * numEqual
        sum = sum + Math.log(omega) * numNotEqual

        if(debugging <= 1) {
            err.println "EpsilonSum: return ${sum}"
        }
        return sum
    } // EpsilonSum
} // EpsilonSumScript
