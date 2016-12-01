FROM ubuntu:latest
MAINTAINER Dave Roe <droe@nmdp.org>

# install stuff
RUN apt-get update && apt-get install -qyy bwa curl git make \
     unzip samtools maven subversion gzip \
     openjdk-7-jdk groovy bowtie2 \
  && cd /opt && git clone "https://github.com/lindenb/jvarkit.git" \
  && cd /opt/jvarkit && make sam2tsv \
  && cd /opt && git clone "https://github.com/google/guava" \
  && cd /opt/guava/guava && mvn install \
  && cd /opt && svn checkout https://svn.code.sf.net/p/dishevelled/code/tags/dsh-commandline-1.1 dsh-commandline \
  && cd /opt/dsh-commandline/ && mvn install \
  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash

# env vars
ENV PATH /opt/jvarkit/dist-1.133:$PATH
ENV PATH /opt/bin:$PATH
ENV CLASSPATH /opt/jvarkit/dist-1.133/sam2tsv.jar:$CLASSPATH
ENV CLASSPATH /opt/guava/guava/target/guava-19.0-SNAPSHOT.jar:$CLASSPATH
ENV CLASSPATH /opt/dsh-commandline/target/dsh-commandline-1.1.jar:$CLASSPATH
ENV NXF_OPTS "-Xms1G -Xmx4G"

# kharsh source
ADD *.nf /opt/kharsh/
ADD tutorial /opt/kharsh/tutorial/
ADD src /opt/kharsh/src/
ENV PATH /opt/kharsh/src:$PATH
ENV CLASSPATH /opt/kharsh/src/kharsh_0.1.jar:/opt/kharsh/src/CsvFileReader_2.1.jar:$CLASSPATH

CMD ["/bin/bash"]
