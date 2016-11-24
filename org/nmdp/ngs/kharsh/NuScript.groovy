package org.nmdp.ngs.kharsh

/*
 * Implements the Nu equation in Equation 3 of the HARSH paper: 
 * the edge-potential function with three outcomes depending on 
 * if the read is the opposite of a predicted haplotype.
 *
 * r is the current assignment of the read relative to haplotype {-1,1}.
 *
 * rj indicates one allele (1) or another (-1) in a read.
 *
 * hi indicates one allele (1) or another (-1) in an estimated haplotype.
 *
 * @author Dave Roe
 * @version $Id: NuScript.groovy 24974 2015-08-17 00:08:37Z droe $
 */

class NuScript {
    static err = System.err
    static int debugging = 5
    static Double Nu(int r, int j, int h, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Nu(j=${j}, h=${h}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if ((r == 1) && (j == h)) { 
            ret = ret + Math.log(epsilon)
        } else if ((r == 1) && (j != h)) {  
            ret = ret + Math.log(1 - epsilon)
        } else { // r == -1
            ret = 0
        }

        if(debugging <= 1) { 
            err.println "Nu: return ${ret}"
        }
        return ret
    } // Nu
} // NuScript
