package org.nmdp.ngs.kharsh

/*
 * Implements the Theta equation in Equation 3 of the HARSH paper.
 * It is an edge-potential function with three outcomes.
 *
 * hi indicates one allele (1) or another (-1) in an estimate haplotype.
 *
 * rj indicates one allele (1) or another (-1) in a read.
 * 
 * e.g., 
 * ln(1 - 0.002) = -0.002; exp(-0.002) = 0.998
 * ln(0.002) = -6.21; exp(-6.21) = 0.002
 *
 * @author
 * @version $Id: ThetaScript.groovy 22961 2015-06-06 01:12:40Z droe $
 */
class ThetaScript {
    static err = System.err
    static int debugging = 6

    static Double Theta(int rj, int hi, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Theta(rj=${rj}, hi=${hi}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if (rj == hi) { 
            ret = ret + Math.log(1 - epsilon)
        } else if (rj != hi) {  
            ret = ret + Math.log(epsilon)
        } else { // rj == -1
            ret = 0
        }

        if(debugging <= 1) { 
            err.println "Theta: return ${ret}"
        }
        return ret
    } // Theta
} // ThetaScript
