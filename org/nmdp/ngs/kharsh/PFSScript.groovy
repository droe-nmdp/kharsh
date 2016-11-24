package org.nmdp.ngs.kharsh

/*
 * Sample reference haplotype on fixed read origins and predicted haplotype.
 * Equation 9 from the HARSH paper.
 *
 * Also, see Figure 15.4 in algorithms.pdf.
 * http://aima.cs.berkeley.edu/
 *
 * All indexes are 0 based.
 *
 * H is a vector representing a predicted haplotype. Each column is one
 * of M indicator variables in {-1,1} where M is the number of variants.
 *
 * S matrix of known/reference haplotypes. Each row is a reference
 * sequence and each column is a variant value {-1,1}. M variants.
 *
 * omega is the probability of haplotype copying 'error'
 * rho is the population recombination rate
 *
 * Return the most likely reference sequence, its reference indexes, and
 * the normalized variant score/probability.
 *
 * @author Dave Roe
 * @version $Id: PFSScript.groovy 25438 2015-08-29 13:03:37Z droe $
 * @todo: modularize this
 */

class PFSScript {
    static err = System.err
    static int debugging = 4 // higher value means more logging
        
    static ArrayList PFS(int[] H, int[][] S, Double omega, Double rho) {
        if(debugging <= 3) {
            err.println "PFS()"
            err.println "H=${H}"
            //err.println "S=${S}"
        }

        int N = S.length // number of reference haplotypes
        int L = H.length // number of variants

        // return sequence
        int[] retS = new int[L]
        // return index to reference haplotype for each snp
        int[] retI = new int[L]

        ArrayList<Double> forward = new ArrayList<Double>(N)
        // calculate the initial state for first variant
        // using just emission/epsilon
        (0..N-1).each { ri -> // each reference haplotype (row index)
            // calculate emission potential
            int siSNP = H[0]
            int refiSNP = S[ri][0]
            if(debugging <= 2) {
                err.println "\nPFS(forward0): H[0] snp 0=${siSNP}, " +
                    "ref ${ri} snp 0=${refiSNP}"
            }
            Double e = Math.exp(EpsilonScript.Epsilon(siSNP, refiSNP, omega))
            forward.add(e)
        } // each reference haplotype

        if(debugging <= 2) {
            err.println "PFS: after initial variant, forward=" + forward
        }
        
        // sample on the first variant
        ProbabilitySampler probabilitySampler = new ProbabilitySampler(forward)
        //err.println "before sample(), forward[0]=" + forward[0]//remove(todo)
        Integer sIdx = probabilitySampler.sample()
        /*
        err.println "after sample(), forward[0]=" + forward[0]//remove(todo)
        err.println "after sample(), prob[0]=" + probabilitySampler.getProbability(0)//remove(todo)
        err.println "randomly sampled ${sIdx}" //todo:remove
        err.println "randomly sampled forward=" + forward[sIdx]//todo:remove
        err.println "randomly sampled prob=" + probabilitySampler.getProbability(sIdx)//todo:remove
        sIdx = 15 // 02:06:02 (todo: testing)
        if(H[0] == 1) {
            sIdx = 107 // 07:01:01:01
        }
        */
        if(debugging <= 2) { 
            err.println "ref hap ${sIdx}=" + S[sIdx]
            err.println "forward ${sIdx}=" + forward[sIdx]
            err.println "prob ${sIdx}=" + probabilitySampler.getProbability(sIdx)
        }

        // forward score is sum of emission scores per ref hap (N of them)
        (1..L-1).each { si -> // each variant index
            int siSNP = H[si]
            (0..N-1).each { ri -> // each reference haplotype (row index)
                // get the previous V value for this haplotype at the previous SNP
                Double Vminus1 = probabilitySampler.getProbability(ri)
                // tau: switch score this reference haplotype relative to the one
                // that was sampled for the previous variant
                Double t = TauScript.Tau(sIdx, ri, N, rho)

                // calculate emission potential
                int refsiSNP = S[ri][si]
                Double e = EpsilonScript.Epsilon(siSNP, refsiSNP, omega)
                Double expSum = Math.exp(t + e)
                Double V = Vminus1 * expSum
                if(debugging <= 2) {
                    err.println "\nPFS(forward1+): ref index ${ri}"
                    err.println "PFS: H[${si}]=${siSNP}, " +
                        "ref snp=${refsiSNP}"
                        err.println String.sprintf(
                            "PFS: Vminus1=%1.5f, exp(%1.5f+%1.5f)=%1.5f",
                            Vminus1, t, e, expSum)
                    err.println "PFS: V=${V}"
                }
                forward.set(ri, V)
             } // each reference haplotype
            // normalize and sample haplotype index
            probabilitySampler = new ProbabilitySampler(forward)
            sIdx = probabilitySampler.sample()
            if(debugging <=2) {
                err.println "PFS: variant ${si} sampled hap index ${sIdx}"
                err.println "PFS: forward=" + forward
            }
        } // each variant
        // end forward
        
        // normalize and sample on the last variant (index L-1)
        probabilitySampler = new ProbabilitySampler(forward)
        sIdx = probabilitySampler.sample()
        Double totalSum = probabilitySampler.getProbability(sIdx)

        retS[L-1] = S[sIdx][L-1]
        retI[L-1]  = sIdx

        if(debugging <= 3) {
            err.println "PFS: after forward"
            err.print "PFS: sampled haplotype index ${sIdx} for last variant: "
            err.println retS[L-1]
            //err.println "haplotypes and scores: " + forward
        }

        // calculate backward starting at variant L-1 (index L-2)
        (L-2..0).each { si -> // each variant (column)
            ArrayList<Double> backward = new ArrayList<Double>(N)
            (0..N-1).each { ri -> // each reference haplotype (row index)
                // get the previous V value for this haplotype at the previous SNP
                Double Vminus1 = probabilitySampler.getProbability(ri)
                // tau: switch score this reference haplotype relative to the one
                // that was sampled for the previous variant
                Double t = TauScript.Tau(retI[si+1], ri, N, rho)
                // epsilon: score this variant at this reference haplotype
                // relative to the variant in the predicted hapltype
                int siSNP = H[si]
                int refiSNP = S[ri][si]
                if(debugging <= 2) {
                    err.println "PFS(backward): H[${ri}] snp ${si}=${siSNP}, " +
                        "ref snp=${refiSNP}"
                }
                Double e = EpsilonScript.Epsilon(siSNP, refiSNP, omega)
                // new score
                backward[ri] = Vminus1 * Math.exp(t + e)
            } // each reference haplotype

            // normalize and sample on the last variant
            probabilitySampler = new ProbabilitySampler(backward)
            sIdx = probabilitySampler.sample()
            totalSum += probabilitySampler.getProbability(sIdx)

            retS[si] = S[sIdx][si]
            retI[si]  = sIdx

            if(debugging <= 2) {
                err.println "PFS: done with backward for index si=${si}"
                err.print "PFS: sampled haplotype index ${sIdx}, variant:"
                err.println retS[si]
                //err.println "haplotypes and scores: " + backward
            }
        } // each variant

        if(debugging <= 3) {
            err.println "PFS: return"
            err.println "PFS: H = " + H
            err.println "PFS: retS = " + retS
            err.println "PFS: retI = " + retI
            err.println "PFS: avg totalSum = " + totalSum/L
        }
        return [retS, retI, totalSum/L]
    } // PFS
} // PFSScript
