package org.nmdp.ngs.kharsh

/*
 * Calculate Epsilon as defined on page 2248 of the HARSH paper.
 *
 * Epsilon is an edge potential representing "'haplotype copying', which
 * is motivated by that the predicted haplotype is a mosaic of reference
 * haplotypes with a small number of differences.
 * 
 * It returns a 0 or negative number, whose exponent is taken
 * as the score.
 *
 * hSNP is the genotype value from the predicted haplotype
 *
 * sSNP is genotype value from the known/reference haplotype
 *
 * omega represents a novelty rate? (the extent to which the predicted
 * haplotype is expected to differ from reference?). A smaller omega
 * gives a bigger difference than a larger; therefore, a smaller omega
 * (e.g., 0.001) will relatively encourage differences, and a 
 * larger one (e.g., 0.1) will relatively discourage differences. (?)
 * 
 * harsh uses a default of 0.002.
 *
 *  omega = 0.1
 *    same = -0.1053605157
 *    different = -2.302585093
 *    difference = 2.1972245
 *
 *  omega = 0.0001
 *    same = -0.0001000050003
 *    different = -9.210340372
 *
 * @author Dave Roe
 * @version $Id: EpsilonScript.groovy 24979 2015-08-17 00:17:56Z droe $
 */

class EpsilonScript {
    static err = System.err
    static int debugging = 2
        
    static Double Epsilon(int hSNP, int sSNP, Double omega) {
        if(debugging <= 1) {
            err.println "Epsilon(hSNP=${hSNP}, sSNP=${sSNP}, omega=${omega})"
        }

        Double e = null
        if((hSNP == null) || (hSNP == 0) || (sSNP == null) || (sSNP == 0)) {
            e = -10  // a very low match
        } else if(hSNP == sSNP) { 
            e = Math.log(1 - omega)
        } else { 
            e = Math.log(omega)
        }

        if(debugging <= 1) {
            err.println "Epsilon: return ${e}"
        }
        return e
    } // Epsilon
} // EpsilonScript
