package org.nmdp.ngs.kharsh

/*
 * Sample haplotype H for fixed reads (R) and reference 
 * haplotype (S): P(H|R,S). It returns the probability that the variant 
 * is from the input haplotype, as opposed to the other haplotype.
 * See equations 7 and 8 in the HARSH paper.
 *
 * To prevent divisions by zero, etc. when taking exponents of negative
 * numbers, implemented via 
 *   1/(1+exp(a-b)) instead of exp(a)/(exp(a)+exp(b)) as written in the 
 * paper.
 * todo: put this back and put into PFR too (and PFS too?)
 *
 * chromosome is 0 or 1, indicating if currentH is the first or second 
 * predicted haplotype.
 *
 * R is a vector of N indicator variables in {-1,1} where N is the number
 * of reads. -1 means is aligns to the first first haplotype, 1 means it
 * aligns to the second haplotype. 
 *
 * currentH is a  vector representing a predicted
 * haplotype. Each column is one of M indicator variables in {-1,1} where
 * M is the number of SNPs.
 *
 * hIndex is the column index into X and S, indicating the sampled variant
 *
 * S is the reference haplotype vector of snps.
 *
 * X is a N by M matrix computed from the read alignment against the reference. 
 * X_ij =0 (unobserved SNP), 1 (major allele), -1 (minor allele). 
 *   Rows are the N reads(j), columns are the N SNPs(i).
 *
 * mu represents the 'heat' of the model (todo: remove)
 * epsilon represents the sequencing error rate; must be > 0.0
 * omega is the haplotype copying 'error'
 *
 * Return the probability that the SNP is from h-bar
 *
 * @author Dave Roe
 * @version $Id: PFHScript.groovy 24975 2015-08-17 00:11:24Z droe $
 */

class PFHScript {
    static err = System.err
    static int debugging = 6

    static Double PFH(int[][] H, int[] hMarkers, int hIndex, int[][] S,
                      int[][] X, int[] R, Double mu, Double epsilon,
                      Double omega) {
        if(debugging <= 1) {
            err.println "PFH(hIndex=${hIndex})"
            err.println "PFH(H=${H})"
        }
        // the alignment column vector for a single SNP
        int otherhIndex = hIndex * -1 + 1
        int[] hAssign = H[hIndex]
        int[] sAssign = S[hIndex]
        int[] hOther = H[otherhIndex]
        int hapMarker = hMarkers[hIndex] // e.g., -1 for H1 and 1 for H2 (reads)
        int otherhapMarker = hMarkers[otherhIndex]
        // get the correct column (variant)
        int[] xi = X.collect { it[hIndex] }
        // get the indexes of the reads assigned to this haplotype
        int[] xHIndexes = [xi, hAssign].transpose().collect { it[0] == it[1] }.findIndexValues { it == true }
        
        // get the indexes of the reads assigned to the other haplotype
        int[] xHOtherIndexes = [xi, hOther].transpose().collect { it[0] == it[1] }.findIndexValues { it == true }

        // the extent to which the variants assigned to this haplotype
        // in X match this haplotype
        Double thetaSumMinus1 = ThetaSumHapsScript.ThetaSumHaps(xHIndexes, xi,
                                                                hAssign,
                                                                epsilon)
        // the extent to which the variants assigned to the other haplotype
        // in X don't match this haplotype
        Double nuSumMinus1 = NuSumHapsScript.NuSumHaps(xHOtherIndexes, xi,
                                                       hAssign, epsilon)
        // the extent to which the prediction matches the reference
        Double epsilon2 = EpsilonScript.Epsilon(hAssign[hIndex],
                                                sAssign[hIndex], omega)

        // numerator
        Double sumMinusTmp = thetaSumMinus1 + nuSumMinus1 + epsilon2
        if(debugging <= 3) {
            err.println "PFH: thetaSumMinus1=${thetaSumMinus1}, nuSumMinus1=${nuSumMinus1}, epsilon2=${epsilon2}"
        }

        // the extent to which the variants assigned to the other haplotype
        // in X match this haplotype
        Double thetaSumPlus1 = ThetaSumHapsScript.ThetaSumHaps(xHOtherIndexes,
                                                               xi, hAssign,
                                                               epsilon)
        // the extent to which the variants assigned to this haplotype
        // in X don't match the other predicted haplotype
        Double nuSumPlus1 = NuSumHapsScript.NuSumHaps(xHIndexes, xi, hOther,
                                                      epsilon)
        Double sumPlusTmp = thetaSumPlus1 + nuSumPlus1 + epsilon2

        if(debugging <= 2) {
            err.println "PFH: thetaSumPlus1=${thetaSumPlus1}, nuSumPlus1=${nuSumPlus1}, epsilon2=${epsilon2}"
        }

        Double difference = sumPlusTmp - sumMinusTmp
        Double differenceExp = Math.exp(difference)

        // denominator
        if(debugging <= 2) { 
            err.println "PFH: sumPlusTmp=${sumPlusTmp}, sumMinusTmp=${sumMinusTmp}"
            err.println "PFH: difference=${difference}, differenceExp=${differenceExp}"
        }

        Double delta = 1/(1+differenceExp);
    
        if(debugging <= 1) { 
            err.println "PFH: return ${delta}, ${sumPlusTmp}, ${sumMinusTmp}"
        }
        Double numerator = Math.exp(sumMinusTmp)
        Double denominator = numerator + Math.exp(sumPlusTmp)
        Double prob =  numerator / denominator

        return prob
    } // PFH
} // PFHScript



