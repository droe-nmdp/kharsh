#!/usr/bin/env groovy

/*
 * Runs KHARSH. 
 * 
 * Input is two matrices (one for the person, one for the reference alleles)
 * and a lot of parameters for the algorithm.
 *
 * Output is IMGT/IPD nomenclature interpretation for both ARS and full gene.
 *
 * e.g., ./kharsh.groovy -m 0.00001 -e 0.01 -o 0.001 -r 0.0001 -c 0.1 -f 800 -g 900 -p input/DRB1-03010101-150201_DRB1_matrix.txt -a input/DRB1-alleles_matrix.txt
 *
 * ./kharsh.groovy -m 1.0E-5 -e 1.0E-4 -o 1.0E-5 -r 0.1 -c 0.1 -f 3000 -g 7500 -p input/2DL4_0080101-00901_8chunks.bwa_2DL4_matrix.txt -a input/2DL4-alleles_matrix2_2.txt 2> output/kharsh_err.txt
 * 
 * e.g., ./harshVCF.groovy -p input/2DL4_0080101-00901_8chunks.bwa_2DL4_matrix.txt -a input/2DL4-alleles_matrix2.txt -o output/2DL4_0080101-00901_8chunks.bwa_2DL4.vcf
 *
 * grep "new best" output/kharsh_err.txt
 *
 * @author Dave Roe
 * @version $Id: kharsh.groovy 24971 2015-08-16 23:58:18Z droe $
 */

import org.dishevelled.commandline.*
import org.dishevelled.commandline.argument.*
import org.nmdp.b12s.kharsh.*

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
out = System.out
err = System.err
IOScript io = new IOScript()

// these are just defaults and shouldn't be used
Double mu = null // "heat" model parameter (e.g., 0.00001)
Double epsilon = null // sequencing error rate (e.g., 0.01)
Double omega = null // haplotype copying 'error' (e.g., 0.001)
Double rho = null // the population recombination rate (e.g., 0.0001)
Double convergeCutoff = null // (e.g., 0.01)
Integer arsStart = null
Integer arsEnd = null
String personMatrixFile
String alleleMatrixFile

(mu, epsilon, omega, rho, convergeCutoff, arsStart, arsEnd, personMatrixFile,
 alleleMatrixFile) = handleArgs(args)

// X is a N by M matrix computed from the read alignment against the reference.
//   X_ij =0 (not matched to the location),
//   Rows are reads, columns are the snps.
int[][] X = null
ArrayList<String> readNames = null
// S matrix of known/reference haplotypes. Each row is a reference
// sequence and each column is a SNP value {-1,1}. M snp values.
int[][] S = null
ArrayList<String> alleleNames = null
List<Integer> commonSnpIndexes = null
(X, readNames, S, alleleNames, commonSnpIndexes) = io.LoadMatrix(personMatrixFile,
                                                                 alleleMatrixFile)
if(debugging <= 2) {
    err.println "X size=${X.size()}"
    //err.println "(from input)X=" + X
    //new FileWriter(new File("X_out.txt")).println X
    err.println "S size=${S.size()}"
    err.println "${S[0].size()} variants"
    err.println "(from input)S=" + S
    //err.println "readNames=${readNames.join(' ')}"
    err.println "alleleNames=${alleleNames.join(' ')}"
}
if(debugging <= 3) {
    err.println "common snps:"
    err.println commonSnpIndexes.join(",")
}

Algorithm2Script alg2 = new Algorithm2Script()
ArrayList<ArrayList> alg2List = alg2.BestOf(S, X, mu, epsilon, omega,
                                            rho, convergeCutoff,
                                            arsStart, arsEnd,
                                            commonSnpIndexes)

ArrayList<ArrayList<Integer>> alg2H = alg2List.get(0)
Float algP = alg2List.get(1)[0]
ArrayList<Integer> predH1 = alg2H.get(0)
ArrayList<Integer> predH2 = alg2H.get(1)

if(debugging <= 3) {
    err.println "prediction:"
    err.println predH1
    err.println predH2
}
// interpret
String interp = alg2.Interpret(S, alleleNames, arsStart, arsEnd, predH1, predH2,
                              commonSnpIndexes)
out.println "ARS: ${interp}"
arsStart = -1; arsEnd = -1
interp = alg2.Interpret(S, alleleNames, arsStart, arsEnd, predH1, predH2,
                        commonSnpIndexes)
out.println "Gene: ${interp}"
if(debugging <= 3) {
    err.println "ARS: ${interp}"
    err.println "Gene: ${interp}"
}
// end main
 
/*
 * handleArgs
 * 
 * @param args the command line arguments
 * @return list of command line arguments in a usable form; use -h for more info
 */
def ArrayList handleArgs(String[] args) {
    USAGE = "kharsh.groovy [args]";
    Switch help = new Switch("h", "help", "display help message");

    muArg = new DoubleArgument("m", "mu", "'heat' model parameter", true)
    epsilonArg = new DoubleArgument("e", "epsilon", "sequencing error rate",
                                    true)
    omegaArg = new DoubleArgument("o", "omega", "haplotype copying 'error'",
                                  true)
    rhoArg = new DoubleArgument("r", "rho",
                                "the population recombination rate", true)
    convergeCutoffArg = new DoubleArgument("c", "convergeCutoff",
                                           "% delta to determine convergence",
                                           true)
    arsStartArg = new IntegerArgument("f", "arsStart",
                                      "index of start of ARS region", true)
    arsEndArg = new IntegerArgument("g", "argsEnd",
                                    "index of end of ARS region", true)
    personFileArg = new FileArgument("p", "person", "person matrix file",
                                     true)
    alleleFileArg = new FileArgument("a", "allele", "allele matrix file",
                                     true)
    ArgumentList arguments =
        new ArgumentList(help, muArg, epsilonArg, omegaArg, rhoArg,
                         convergeCutoffArg, arsStartArg, arsEndArg,
                         personFileArg, alleleFileArg)
    CommandLine commandLine = new CommandLine(args);
    try {
        CommandLineParser.parse(commandLine, arguments);
        if (help.wasFound()) {
            Usage.usage(USAGE, null, commandLine, arguments, System.out);
            System.exit(-2)
        }
    } catch (CommandLineParseException e) {
        Usage.usage(USAGE, e, commandLine, arguments, System.err);
        System.exit(-1)
    }

    return [muArg.getValue(), epsilonArg.getValue(), omegaArg.getValue(),
            rhoArg.getValue(), convergeCutoffArg.getValue(),
            arsStartArg.getValue(), arsEndArg.getValue(),
            personFileArg.getValue(), alleleFileArg.getValue()]
} // handleArgs