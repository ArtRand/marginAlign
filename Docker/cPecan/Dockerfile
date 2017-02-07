# Basics
FROM gcc:5.4

MAINTAINER Art Rand, arand@soe.ucsc.edu

# apt-get installs
RUN apt-get update && apt-get install -y git make zlib1g-dev g++
WORKDIR /home/

# sonLib
RUN git clone https://github.com/benedictpaten/sonLib.git
WORKDIR /home/sonLib/
RUN make

# cPecan
WORKDIR /home/
RUN git clone https://github.com/benedictpaten/cPecan.git
WORKDIR /home/cPecan/
RUN make


# prepare
#RUN ls -al /home/sonLib/bin
RUN mkdir /opt/cPecanRealign && \
          cp  /home/sonLib/bin/cPecanRealign /opt/cPecanRealign/

COPY runner.py /opt/cPecanRealign/runner
RUN chmod +x /opt/cPecanRealign/runner
COPY wrapper.sh /opt/cPecanRealign

RUN mkdir /data
WORKDIR /data

# set PATH variables
ENV PATH "/opt/cPecanRealign/:$PATH"

ENTRYPOINT ["sh", "/opt/cPecanRealign/wrapper.sh"]
