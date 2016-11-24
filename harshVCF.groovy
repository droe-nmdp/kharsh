#!/usr/bin/env groovy

/*
 * Makes a vcf of the alleles from the reference alleles on _the_ 
 * reference allele.
 * 
 * It requires the kharsh jar and Dishevelled CommandLine:
 *   http://www.dishevelled.org/commandline/
 *
 * Input is the matrix file (e.g., test1-alleles_2DL4_matrix_1_uniq_3.txt)
 *
 * Output is a vcf with chrom#, snpID, position, position, referenceAllele, altAllele
 * From converter.py:
 *       data.chrom + "\tSNP" + str(count) + "\t" + 
 *			   str(data.pos) + "\t" + str(data.pos) +  "\t" + 
 *			   data.ref + "\t" + data.alt + "\n"
 *
 * http://genetics.cs.ucla.edu/harsh/manual.html
 *
 * e.g., harshVCF.groovy -p input/2DL4_0080101-00901_8chunks.bwa_2DL4_matrix.txt -a input/2DL4-alleles_matrix2.txt -i input/2DL4_0080101-00901_8chunks.bwa_2DL4_variants.txt -o output/2DL4_0080101-00901_8chunks.bwa_2DL4.vcf
 *
 * @author Dave Roe
 * @version $Id: harshVCF.groovy 32774 2016-05-02 20:18:39Z droe $
 */

import org.dishevelled.commandline.*
import org.dishevelled.commandline.argument.*
import org.nmdp.b12s.kharsh.*

// things that may change per run
debugging = 4 // TRACE=1, DEBUG=2, INFO=3
out = System.out
err = System.err

// cast to string
String personFile
String alleleFile
FileReader vInFile
String outFile
(personFile, alleleFile, vInFile, outFile) = handleArgs(args)
IOScript io = new IOScript()
io.makeHarshVcf(personFile, alleleFile, vInFile, outFile)

/*
 * handleArgs
 * 
 * @param args the command line arguments
 * @return list of command line arguments in a usable form; use -h for more info
 */
def ArrayList handleArgs(String[] args) {
    USAGE = "kharsh.groovy [args]";
    Switch help = new Switch("h", "help", "display help message");

    personFileArg = new FileArgument("p", "person", "person matrix file",
                                     true)
    alleleFileArg = new FileArgument("a", "allele", "allele matrix file",
                                     true)
    varInFileArg = new FileArgument("i", "variants input",
                                    "variants input file",
                                    false)
    outFileArg = new FileArgument("o", "out", "output VCF file",
                                  true)

    ArgumentList arguments = new ArgumentList(help, personFileArg, alleleFileArg,
                                              varInFileArg, outFileArg)
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

    return [personFileArg.getValue(), alleleFileArg.getValue(),
            varInFileArg.getValue(), outFileArg.getValue()]
} // handleArgs

