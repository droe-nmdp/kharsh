package org.nmdp.ngs.kharsh

/*
 * Implements the Theta equation in Equation 3 of the HARSH paper: 
 * the edge-potential function with three outcomes depending on 
 * if the read is mapped to the predicted haplotype.
 * 
 * This assumes the reads have been assigned to the haplotype.
 *
 * x indicates an allele value in a read
 *
 * hi indicates an allele value in a haplotype.
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

    static Double Theta(int x, int h, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Theta(x=${x}, h=${h}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if (x == h) { 
            ret = ret + Math.log(1 - epsilon)
        } else {  
            ret = ret + Math.log(epsilon)
        }

        if(debugging <= 1) { 
            err.println "Theta: return ${ret}"
        }
        return ret
    } // Theta
} // ThetaScript
