package org.nmdp.ngs.kharsh

/*
 * Test equality of two numbers +/- a cutoff.
 *
 * v1 a number; the new values
 * v2 a number; the old value
 * a cutoff value
 *
 * Returns 1 if the v1 >= v2 and the difference is less than
 * or equal to the percentage in cutoff.
 *
 * @author Dave Roe
 * @version $Id: ConvergedScript.groovy 22961 2015-06-06 01:12:40Z droe $
 */

class ConvergedScript {
    static err = System.err
    static int debugging = 3

    static Boolean Converged(Double v1, Double v2, Double cutoff) {
        if(debugging <= 1) { 
            err.println "Converged(v1=${v1}, v2=${v2}, cutoff=${cutoff})"
        }
        Boolean bool = false // assume false
        if((v1 != 0) || (v2 != 0)) { 
            Double fraction = Math.abs(v1-v2)/v1
            if(debugging <= 2) { 
                err.println "fraction=${fraction}, cutoff=${cutoff}"
            }
            if((v1 >= v2) && (fraction <= cutoff)) { 
                bool = true
            }
        } else { 
            bool = true
        }
        if(debugging <= 1) { 
            err.println "Converged: return ${bool}"
        }
        return bool
    } // Converged
} // ConvergedScript
