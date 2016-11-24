package org.nmdp.ngs.kharsh

/*
 * Calculate Tau as defined on page 2247-2248 of the HARSH paper. 
 * It starts as step 2 above Equation 6.
 *
 * Tau "model(s) the transition probability in (the) haplotype copying
 * model (Li and Stephens, 2003)".
 * 
 * s1 and s2 are indexes-into or indicators-of two reference haplotypes;
 * i.e., if s1==s2 they are the same haplotype; if s1!=s2 they
 * are different haplotypes
 * 
 * N is # of reference haplotypes
 *
 * rho is the population recombination rate
 * A smaller rho returns a larger tau and therefore deemphasizes recombination.
 * e.g., 
 *   rho 0.1, N=10
 *   same = 0.9910448503
 *   different = 0.00099501663
 *   
 *   rho 0.0001
 *   same = 0.999991
 *   different = 0.000001
 *   
 *   rho = 0.9
 *   same = 0.9225380668
 *   different 0.00860688147
 *
 * @author Dave Roe
 * @version $Id: TauScript.groovy 25293 2015-08-22 14:19:17Z droe $
 */

class TauScript {
    static err = System.err
    static int debugging = 2
        
    static Double Tau(int s1, int s2, int N, Double rho) {
        if(debugging <= 1) {
            err.println "Tau: s1=${s1}, s2=${s2}, N=${N}, rho=${rho}"
        }
        Double val = 0
        if(s1 == s2) { 
            val = Math.exp(-rho/N) + (1-Math.exp(-rho/N))/N
        } else { 
            val = (1-Math.exp(-rho/N))/N
        }

        if(debugging <= 1) {
            err.println "Tau: return ${val}"
        }
        return val
    } // Tau
} // TauScript
