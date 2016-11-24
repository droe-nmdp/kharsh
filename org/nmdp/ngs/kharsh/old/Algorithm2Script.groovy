package org.nmdp.ngs.kharsh

/*
 * Use Gibbs sampling of reads, two predicted haplotypes, and a
 * collection of reference haplotypes to predict a pair of chromosome 
 * given a vector of reads mapped to a reference haplotype, a matrix
 * of reads mapped to a collection of reference haplotypes, and a 
 * matrix of SNPs and reads.
 *
 * Sample until both read and haplotype probabilities converge
 * to within convergeCutoff or a certain number of rounds.
 *
 * @todo get rid of the converge cutoff stuff? doesn't really do anything
 * right now
 *
 * @author Dave Roe
 * @version $Id: Algorithm2Script.groovy 25302 2015-08-23 16:13:30Z droe $
 *
 */

class Algorithm2Script {
    // the higher, the less output
    static Integer debugging = 1
    static Integer maxSamples = 1 // # of calls to Algorithm2()
    static Integer maxRounds = 1 // # of iterations of alg2 per sample
    static Random rand = new Random()
    static err = System.err

    ArrayList<ArrayList> Algorithm2(int[][] S, int[][] X, Double mu,
                                    Double epsilon, Double omega, Double rho,
                                    Double convergeCutoff) {
        Integer N = X.size()   // number of reads
        Integer M = X[0].size()   // number of SNPs
        // Z is a normalization constant to ensure the sum of P() for all R, H
        // equal 1
        Double Z = 1 //TODO
        Integer burnInRounds = 1 // todo: use a coverage calculation here(?)    
        // probabilites for each read
        int[] rProbs = new int[N]

        // Step 1: Randomly initialize haplotype H to SNPs 1 or -1
        // and randomly initialize read assignment to 1 or -1.
        int[][] H
        int[][] returnH
        int[] R
        (H, R) = Algorithm2Init(X, N, M)

        // randomly initialize the two reference haplotypes; contents
        // of S1 and S2 are varient values; contents of SI1 and SI2 are
        // indexes into a reference haplotype
        ArrayList sampledS, sampledSI
        (sampledS, sampledSI) = Algorithm2InitRefs(S, M)

        def H1 = H[0];             def H2 = H[1]
        def S1 = sampledS.get(0);      def S2 = sampledS.get(1)
        def SI1 = sampledSI.get(0);    def SI2 = sampledSI.get(1)
        double[] sampledScores = new double[2]

        int i = 0 // index for number of rounds
        boolean converged = false
        Double pPrevious = 0 // overall probability from the previous round
        Double p = PScript.P(R, H1, H2, S1, S2, SI1, SI2, S.size(), X, mu, 
                             epsilon, omega, rho, Z)
        Double returnP=p
        if(debugging <= 3) {
            err.println String.sprintf("Algorithm2: after init, p=%1.5f\n", p)
        }
        ArrayList pList = new ArrayList() // history of the overall probability
        pList.add(p)
        returnH = H
        // sample until both reads and haplotypes converge or maxRounds is reached
        while((i < burnInRounds) || ((i < maxRounds) && (converged == false))) {
            H1 = H[0]
            H2 = H[1]
            if(debugging <= 3) { 
                err.println String.sprintf("\nstarting round %d, p=%1.5f\n", i, p)
            }
            // Step 2: For fixed haplotype H, sample read origin R
            // see equation 6
            (0..N-1).each { rIndex ->        // each of the N reads
                Double prob = PFRScript.PFR(H, X, rIndex, mu, epsilon)
                Double r = (Double)Math.random()
                if(debugging <= 2) {
                    err.println "after PFR: rIndex=${rIndex}, prob=${prob}, " +
                        "r=${r}"
                }
                // assuming eq 6 should be P(rj) = -1
                //todo: is this right?
                if(prob > r) { 
                    R[rIndex] = -1
                } else { 
                    R[rIndex] = 1
                }
                rProbs[rIndex] = prob
            } // each read

            (0..1).each { chromosome ->
                int[] currentH = H[chromosome]
                // Step 3: For fixed haplotype H, sample haplotype reference S
                int[] tmpS; int[] tmpSI; double tmpSum
                (tmpS, tmpSI, tmpSum) = PFSScript.PFS(currentH, S, omega, rho)
                sampledS.set(chromosome, tmpS)
                sampledSI.set(chromosome, tmpSI)
                sampledScores[chromosome] = tmpSum
            } // PFS for each chromosome

            if(debugging <= 3) {
                err.println "sampledS=" + sampledS
                err.println "sampledSI=" + sampledSI
                err.println "sampledScores=" + sampledScores
            }

            // Step 4: For fixed read origin R and haplotype reference S,
            // sample haplotype H; see equations 7&8
            (0..M-1).each { hIndex -> // each snp
                Double h1Prob = 0
                Double h2Prob = 0
                (0..1).each { chromosome -> 
                    int[] currentH = H[chromosome]
                    def currentS = sampledS.get(chromosome)
                    def (prob, sumPlus1, sumMinus1) =
                        PFHScript.PFH(chromosome, R, currentH, X, hIndex,
                                      currentS, mu, epsilon, omega)
                    // sample on the returned probability
                    Double r = (Double)Math.random()
                    if(debugging <= 2) {
                        err.println "after PFH: hIndex=${hIndex}, chromosome=${chromosome}, " +
                            "prob=${prob}, r=${r}"
                    }
                    if(prob > r) {
                        currentH[hIndex] = -1
                    } else { 
                        currentH[hIndex] = 1
                    }
                } // each chromosome
            } // each SNP

            H[0] = H1
            H[1] = H2
            if(debugging <= 2) {
                err.println "H[0]=" + H[0]
                err.println "H[1]=" + H[1]
            }
            //todo(remove)
//todo(remove)            if((H[0] == [1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1]) || (H[1] == [1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1])) {//todo(testing)
//todo(remove)                err.println "*match*"
//todo(remove)                if((H[0] == [-1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1]) || (H[1] == [-1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1])) {
//todo(remove)                    err.println "**double match**"
//todo(remove)                }
                //System.exit(1)
//todo(remove)            }
            p = PScript.P(R, H1, H2, S1, S2, SI1, SI2, S.size(), X, mu, 
                          epsilon, omega, rho, Z)
            if(debugging <= 2) {
                err.println "p=${p}"
            }

            converged = ConvergedScript.Converged(p, pPrevious, convergeCutoff)
            if(debugging <= 2) { 
                err.println "pList max=" + pList.max()
            }
            if(p > pList.max()) {
                returnH = H
                returnP = p
            }
            pList << p
            pPrevious = p
            i = i + 1
            if(debugging <= 3) {
                err.println String.sprintf("\nend round %d, p=%1.5f\n", i, p)
            }
            //System.exit(1)//todo
        } // each round

        if(debugging <= 5) {
            if(i == maxRounds) { 
                err.println String.sprintf("\ndone: reached max %d rounds; p=%1.5f", i, returnP)
            } else { 
                err.println String.sprintf("\ndone: reached %d rounds; p=%1.5f", i, returnP)
            }
            err.println "Algorithm2: returning " + returnH[0]
        }
        pList = new ArrayList() // history of the overall probability
        pList.add(returnP)
        return [returnH, pList]
    } // Algorithm2

    /*
     * Algorithm2Init
     *
     * Randomly initialize haplotype H to snps 1 or -1. i.e., one
     * allele or the other.
     * Randomly initialize read assignment R to 1 or -1. i.e., one
     * chromosome or the other.
     *
     * @param X matrix mapping reads to predicted haplotypes;
     *          should be limited to common heterozygous snps 
     *          in the individual
     * @param N number of reads (# of rows in X)
     * @param M number of SNPs (# of columns in X)
     * @return ArrayList of randomly predicted haplotypes (H) 
     *         and read assignments (R).
     */
    ArrayList<int[][]> Algorithm2Init(int[][] X, int N, int M) {
        int[] R = new int[N]
        (1..N).each { i ->
            R[i-1] = (int)rand.nextInt(2).intValue()
        }
        def zeroIndexes = R.findIndexValues { it == 0 }
        zeroIndexes.each { R[it] = -1 }
        int[] H1 = new int[M]
        (1..M).each { i ->
            H1[i-1] = (int)rand.nextInt(2).intValue()
        }
        zeroIndexes = H1.findIndexValues { it == 0 }
        zeroIndexes.each { H1[it] = -1 }
        int[] H2 = H1.collect{ it * -1 }
        int[][] H = [[], []]
        H[0] = H1
        H[1] = H2
        /*todo (remove) 
        H[0] = [-1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1] // C*02:06:02//todo
        H[1] = [1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1] // C*07:01:01//todo
        */
        if(debugging <= 3) {
            err.println "\ninitialized chromosomes:"
            err.println H
        }
        if(debugging <= 2) {
            err.println "\ninitialized R:"
            err.println R
        }

        return [H, R]
    } // Algorithm2Init

    /*
     * Algorithm2InitRefs
     *
     * Randomly pick two reference haplotypes from the set.
     *
     * @param S a matrix of reference haplotypes, one per row
     * @param M the number of SNPs per reference haploytpe (# of columns)
     * @return two matrices 
     *         1) two reference haplotypes, rows=haplotypes, columns=variants
     *         2) two index sets, rows=haplotypes, 
     *                            columns=indexes to ref hap for each variant
     */
    ArrayList<int[][]> Algorithm2InitRefs(int[][] S, int M) {
        // hap1
        def numRefHaps = S.size()
        def index = (int)rand.nextInt(numRefHaps)
        def S1index = index
        def sampledS1 = S[S1index]
        // all variants index to S1index ref hap
        def sampledS1I = (new int[M]).collect { it = S1index }

        // hap2
        while((index = (int)rand.nextInt(numRefHaps)) == S1index) {
            ; // no-op
        }
        def S2index = index
        def sampledS2 = S[S2index]
        // all variants index to S2index ref hap
        def sampledS2I = (new int[M]).collect { it = S2index }

        def sampledS = [][]
        def sampledSI = [][]
        sampledS[0] = sampledS1
        sampledS[1] = sampledS2
        sampledSI[0] = sampledS1I
        sampledSI[1] = sampledS2I

        return [sampledS, sampledSI]
    } // Algorithm2InitRefs

    /*
     * BestOf
     *
     * Run Algorithm2 a certain number of times, and return the answer
     * with the highest probability.
     *
     * The input and output is the same as Algorithm2().
     */
    ArrayList<ArrayList> BestOf(int[][] S, int[][] X, Double mu, Double epsilon,
                                Double omega, Double rho, Double convergeCutoff,
                                Integer arsStart, Integer arsEnd,
                                ArrayList<Integer> commonSnpIndexes) {
        if(debugging <= 1) {
            err.println "BestOf()"
        }
        ArrayList<ArrayList> alg2List = Algorithm2(S, X, mu, epsilon, omega,
                                                   rho, convergeCutoff)
        ArrayList<ArrayList<Integer>> alg2H = alg2List.get(0)
        Float algP = alg2List.get(1)[-1]
        ArrayList bestAlg2List = alg2List
        ArrayList<Double> pList = alg2List.get(1)
        Double bestProb = pList[-1] // from the first round
        if(debugging <= 4) { 
            err.println "*initial best *${bestProb}"
        }
        
        ArrayList<Integer> H1 = ShortenToARS(Arrays.asList(alg2H.get(0)),
                                             arsStart, arsEnd,
                                             commonSnpIndexes)
        ArrayList<Integer> H2 = ShortenToARS(Arrays.asList(alg2H.get(1)),
                                             arsStart, arsEnd,
                                             commonSnpIndexes)
        int c = 1 // the first run was passed in
        while(c < maxSamples) {
            c++;
            //Thread.sleep(4000) // debugging
            alg2List = Algorithm2(S, X, mu, epsilon, omega, rho, convergeCutoff)
            alg2H = alg2List.get(0)
            algP = alg2List.get(1)[-1]
            if(algP > bestProb) {
                if(debugging <= 4) { 
                    err.println "*new best *${algP} (sample ${c})"
                }
                bestProb = algP
                bestAlg2List = alg2List
            }
                
            H1 = ShortenToARS(Arrays.asList(alg2H.get(0)), arsStart, arsEnd,
                              commonSnpIndexes)
            H2 = ShortenToARS(Arrays.asList(alg2H.get(1)), arsStart, arsEnd,
                              commonSnpIndexes)
            if(debugging <= 3) { 
                err.println "BestOf: c=${c}"
                err.println "BestOf: ARS H1=${H1}"
            }
        } // while less than max rounds

        if(c == maxSamples) {
            alg2List = bestAlg2List
        }
        if(debugging <= 4) { 
            err.println "BestOf: overall best p=${bestProb}"
        }
        if(debugging <= 3) {
            err.println "BestOf: return; ${c} iterations, best=${bestProb}"
        }
        return alg2List
    } // BestOf
    
    /*
     * DupRun
     *
     * Run Algorithm2 until you get the same answer twice.
     *
     * It assumes Algorithm2 returns two haplotypes that are inverses.
     *
     * @todo document input/output
     */
    ArrayList<ArrayList> DupRun(int[][] S, int[][] X, Double mu, Double epsilon,
                                Double omega, Double rho, Double convergeCutoff,
                                Integer arsStart, Integer arsEnd,
                                ArrayList<Integer> commonSnpIndexes) {
        if(debugging <= 1) {
            err.println "DupRun()"
        }
        ArrayList<ArrayList<Integer>> answers = new ArrayList()
        ArrayList<ArrayList> alg2List = Algorithm2(S, X, mu, epsilon, omega,
                                                   rho, convergeCutoff)
        ArrayList<ArrayList<Integer>> alg2H = alg2List.get(0)
        Float algP = alg2List.get(1)[-1]
        ArrayList bestAlg2List = alg2List
        ArrayList<Double> pList = alg2List.get(1)
        Double bestProb = pList[-1]
        
        ArrayList<Integer> H1 = alg2H.get(0)
        ArrayList<Integer> H2 = alg2H.get(1)
        H1 = ShortenToARS(Arrays.asList(H1), arsStart, arsEnd, commonSnpIndexes)
        H2 = ShortenToARS(Arrays.asList(H2), arsStart, arsEnd, commonSnpIndexes)
        int c = 1
        while(!answers.contains(H1) && (c < maxRounds)) {
            // this assumes reciprocity of the two chromosomes
            answers.add(H1)
            answers.add(H2)
            c++;
            //Thread.sleep(4000) // debugging
            alg2List = Algorithm2(S, X, mu, epsilon, omega, rho, convergeCutoff)
            alg2H = alg2List.get(0)
            algP = alg2List.get(1)[-1]

            H1 = Arrays.asList(alg2H.get(0))
            H2 = Arrays.asList(alg2H.get(1))
            H1 = ShortenToARS(Arrays.asList(H1), arsStart, arsEnd, commonSnpIndexes)
            H2 = ShortenToARS(Arrays.asList(H2), arsStart, arsEnd, commonSnpIndexes)
            if(debugging <= 3) { 
                err.println "DupRun: ${answers.size()} answers, c=${c}"
                //err.println "DupRun: H1=${H1}"
            }
        } // while we haven't seen this answer yet

        if(debugging <= 3) {
            err.println "DupRun: return; ${c} iterations, max=${maxRounds}"
        }
        return alg2List
    } // DupRun

    /*
     * Interpret
     *
     * @param startIndex a 0-based start index of variants to include in interpretation
     * @param endIndex a 0-based, end-inclusive end index of variants to
     *                 include in interpretation
     * @todo document
     */
    String Interpret(int[][] S, ArrayList<String> alleleNames,
                     Integer startIndex, Integer endIndex,
                     ArrayList<Integer> H1, ArrayList<Integer> H2,
                     List<Integer> commonSnpIndexes) {
        if(debugging <= 1) {
            err.println "Interpret(startIndex=${startIndex}, endIndex=${endIndex})"
        }
        if((startIndex == -1) || (endIndex == -1)) {
            startIndex = commonSnpIndexes[0]
            endIndex = commonSnpIndexes[-1]
        }
        
        Map<String, Boolean> h1NameMap = [:]
        Map<String, Boolean> h2NameMap = [:]
        (0..S.size()-1).each { key ->
            h1NameMap[alleleNames[key]] = true
            h2NameMap[alleleNames[key]] = true
        }
        Map<String, Boolean> arsMap = [:]
        S.eachWithIndex { s, sIndex -> // s is int[]
            if(sIndex == 0) { 
                err.println "${s.size()} total variants"
            }
            s.eachWithIndex { v, vIndex ->
                if((commonSnpIndexes[vIndex] >= startIndex) &&
                   (commonSnpIndexes[vIndex] <= endIndex)) {
                    arsMap[vIndex] = true
                    if(v != H1[vIndex]) {
                        if(debugging <= 2) {
                            err.println "Interpret: for H1, removing variant ${v}(${vIndex}) for allele index ${sIndex}"
                        }
                        h1NameMap.remove(alleleNames[sIndex])
                    }
                    if(v != H2[vIndex]) {
                        if(debugging <= 2) {
                            err.println "Interpret: for H2, removing variant ${v}(${vIndex}) for allele index ${sIndex}"
                        }
                        h2NameMap.remove(alleleNames[sIndex])
                    }
                } // if the variant is within range
            } // each variant in the list of reference haplotypes
        } // each reference haplotypes

        String h1Str = "NEW"
        String h2Str = "NEW"
        if(h1NameMap.keySet().size() > 0) {
            h1Str = h1NameMap.keySet().sort().join('/')
        }
        if(h2NameMap.keySet().size() > 0) {
            h2Str = h2NameMap.keySet().sort().join('/')
        }
        String ret = "${h1Str}+${h2Str}"
        if(debugging <= 3) {
            err.println "${arsMap.keySet().size()} ARS variants"
        }
        if(debugging <= 1) {
            err.println "Interpret: return ${ret}"
        }
        return ret
    } // Interpret

    ArrayList<Integer> ShortenToARS(List<Integer> h, Integer startIndex,
                                    Integer endIndex,
                                    ArrayList<Integer> commonSnpIndexes) {
        int retIndex = 0
        if(endIndex <= startIndex) {
            return h
        }
        ArrayList<Integer> ret = new ArrayList()
        h.eachWithIndex { v, vIndex ->
            if((commonSnpIndexes[vIndex] >= startIndex) &&
               (commonSnpIndexes[vIndex] <= endIndex)) {
                ret[retIndex++] = v
            }
        } // each variant
        return ret
    } // ShortenToARS
} // Algorithm2Script

