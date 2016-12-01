package org.nmdp.ngs.kharsh

/*
 * Implements the Nu equation in Equation 3 of the HARSH paper: 
 * the edge-potential function with three outcomes depending on 
 * if the read is the opposite of a predicted haplotype.
 *
 * This assumes the reads have been assigned to the *other* haplotype.
 *
 * x indicates an allele value in a read
 *
 * hi indicates an allele value in a haplotype.
 *
 * @author Dave Roe
 * @version $Id: NuScript.groovy 24974 2015-08-17 00:08:37Z droe $
 */

class NuScript {
    static err = System.err
    static int debugging = 5

    static Double Nu(int x, int h, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Nu(x=${x}, h=${h}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if (x == h) { 
            ret = ret + Math.log(epsilon)
        } else {  
            ret = ret + Math.log(1 - epsilon)
        }

        if(debugging <= 1) { 
            err.println "Nu: return ${ret}"
        }
        return ret
    } // Nu
} // NuScript
