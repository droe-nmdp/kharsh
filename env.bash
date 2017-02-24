#!/usr/bin/env bash
# e.g., jar -cvf kharsh_0.1.jar org; cp kharsh.groovy kharsh_0.1.jar ~/git/kharsh/src/
# e.g., export CLASSPATH=$HOME/git/kharsh:$CLASSPATH

export PATH=$PWD/src:$PATH
export PATH=$PWD/src/jvarkit/:$PATH
export CLASSPATH=$PWD/src/jvarkit/sam2tsv.jar:$CLASSPATH
export CLASSPATH=$PWD/src/kharsh_0.1.jar:$CLASSPATH
export CLASSPATH=$PWD/src/CsvFileReader_2.1.jar:$CLASSPATH
export CLASSPATH=$PWD/src/dsh-commandline-1.1.jar:$CLASSPATH
export CLASSPATH=$PWD/src/guava-11.0.2.jar:$CLASSPATH
export PATH=$PWD/src/readsMatrix.groovy:$PATH
export NXF_OPTS="-Xms1G -Xmx4G"
