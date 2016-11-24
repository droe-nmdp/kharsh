package org.nmdp.ngs.kharsh

/*
 * Gibbs sampling for haplotype origin -1 for fixed reads and reference
 * haplotypes.
 * e.g., P(hi = -1 | R, S).
 * Equations 7&8 in the HARSH paper. 
 *
 * To prevent divisions by zero, etc. when taking exponents of negative
 * numbers, implemented via 
 *   1/(1+exp(a-b)) instead of exp(a)/(exp(a)+exp(b)) as written in the 
 * paper.
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
        
    static ArrayList PFH(int chromosome, int[] R, int[] currentH, int[][] X,
                         int hIndex, int[] S, Double mu, Double epsilon,
                         Double omega) {
        if(debugging <= 1) {
            err.println "PFH(hIndex=${hIndex}"
            err.println "PFH(currentH=${currentH})"
            err.println "PFH(S=${S})"
        }
        // the alignment column vector for a single SNP
        int[] xi = X.collect { it[hIndex] }
        int sSNP = S[hIndex]

        int rjMultiplier = (int)-1; // for chromosome 1
        if(chromosome == 1) { //todo?
            rjMultiplier = (int)1
        }

        int snp = (int)-1
        int[][] emptyX = [][]
        //todo(changed from snp to sSNP)
        Double thetaSumMinus1 = ThetaSumHapsScript.ThetaSumHaps(R, xi, emptyX, sSNP.intValue(),
                                                                rjMultiplier.intValue(),
                                                                epsilon)
        //todo(changed from snp to sSNP)
        Double nuSumMinus1 = NuSumHapsScript.NuSumHaps(R, xi, emptyX, sSNP,
                                                       rjMultiplier, epsilon)
        Double epsilon2 = EpsilonScript.Epsilon(snp, sSNP, omega)

        // numerator
        Double sumMinusTmp = thetaSumMinus1 + nuSumMinus1 + epsilon2
//        Double sumMinus1 = Math.exp(sumMinusTmp)
        if(debugging <= 3) {
            err.println "PFH: thetaSumMinus1=${thetaSumMinus1}, nuSumMinus1=${nuSumMinus1}, epsilon2=${epsilon2}"
        }

        snp = (int)1
        //todo(changed from snp to sSNP)
        Double thetaSumPlus1 = ThetaSumHapsScript.ThetaSumHaps(R, xi, emptyX, sSNP,
                                                               rjMultiplier,
                                                               epsilon)
        //todo(changed from snp to sSNP)
        Double nuSumPlus1 = NuSumHapsScript.NuSumHaps(R, xi, emptyX, sSNP,
                                                      rjMultiplier, epsilon)
        epsilon2 = EpsilonScript.Epsilon(snp, sSNP, omega)
        Double sumPlusTmp = thetaSumPlus1 + nuSumPlus1 + epsilon2
        //Double sumPlus1 = Math.exp(sumPlusTmp)

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

        return [delta, sumPlusTmp, sumMinusTmp]
    } // PFH
} // PFHScript
