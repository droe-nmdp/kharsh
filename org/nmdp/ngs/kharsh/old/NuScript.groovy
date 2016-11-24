package org.nmdp.ngs.kharsh

/*
 * Implements the Nu equation in Equation 3 of the HARSH paper.
 * It is an edge-potential function with three outcomes.
 *
 * hi indicates one allele (1) or another (-1) in an estimate haplotype.
 *
 * rj indicates one allele (1) or another (-1) in a read.
 *
 * @author Dave Roe
 * @version $Id: NuScript.groovy 24974 2015-08-17 00:08:37Z droe $
 */

class NuScript {
    static Double Nu(int rj, int hi, Double epsilon) {
        if(debugging <= 1) { 
            err.println "Nu(rj=${rj}, hi=${hi}, epsilon=${epsilon}"
        }
        Double ret = 0.0
        if ((rj == 1) && (rj == hi)) { 
            ret = ret + Math.log(1 - epsilon)
        } else if (rj == 1) {  
            ret = ret + Math.log(epsilon)
        } else { // rj == -1
            ret = 0
        }

        if(debugging <= 1) { 
            err.println "Nu: return ${ret}"
        }
        return ret
    } // Nu
} // NuScript
