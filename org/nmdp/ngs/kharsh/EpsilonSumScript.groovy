package org.nmdp.ngs.kharsh

/*
 * Calculate the sum of epsilons given a predicted haplotype and a
 * reference haplotype. See page 2248 of the HARSH paper.
 *
 * Epsilon is an edge potential representing 'haplotype copying', which
 * is motivated by the idea that the predicted haplotype is a mosaic of 
 * reference haplotypes with a small number of differences. It returns
 * a large negative number (>0) when S & H greatly mismatch, and
 * a small negative number (<0) when S & H greatly match.
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
    static int debugging = 1
        
    static Double EpsilonSum(int[] H, int[] S, Double omega) { 
        Double sum = 0.0
        if(debugging <= 1) {
            err.println "EpsilonSum(omega=${omega})"
        }

        if(debugging <= 2) {
            err.println "EpsilonSum: H=${H}"
            err.println "EpsilonSum: S=${S}"
        }
        
        // compare the variants between H and S
        // ignore the 0 variant values in S
        int numEqual = 0
        [H, S].transpose().collect { 
            if(it[0] == it[1]) {
                numEqual++
            }
        }
        int numSZeros = S.findIndexValues{it == 0}.size()
        int numNotEqual = H.length - numEqual - numSZeros

        if(debugging <= 2) {
            err.println "EpsilonSum: numEqual=${numEqual}, " +
                "numNotEqual=${numNotEqual}, numSZeros=${numSZeros}"
        }

        sum = sum + Math.log(1-omega) * numEqual
        sum = sum + Math.log(omega) * numNotEqual

        if(debugging <= 1) {
            err.println "EpsilonSum: return ${sum}"
        }
        return sum
    } // EpsilonSum
} // EpsilonSumScript
