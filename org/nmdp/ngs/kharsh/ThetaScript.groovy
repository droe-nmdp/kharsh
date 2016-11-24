package org.nmdp.ngs.kharsh

/*
 * Implements the Theta equation in Equation 3 of the HARSH paper: 
 * the edge-potential function with three outcomes depending on 
 * if the read is mapped to the predicted haplotype haplotype indicated
 * by 1 and matches, mapped to haplotype 1 but mismatches, and is mapped
 * to the other predicted haplotype.
 *
 * r is the current assignment of the read relative to haplotype {-1,1}.
 *
 * rj indicates one allele (1) or another (-1) in a read.
 *
 * hi indicates one allele (1) or another (-1) in an estimated haplotype.
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

    static Double Theta(int r, int j, int h, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Theta(r=${r}, h=${h}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if ((r == 1) && (j == h)) { 
            ret = ret + Math.log(1 - epsilon)
        } else if ((r == 1) && (j != h)) {  
            ret = ret + Math.log(epsilon)
        } else { // r == -1
            ret = 0
        }

        if(debugging <= 1) { 
            err.println "Theta: return ${ret}"
        }
        return ret
    } // Theta
} // ThetaScript
