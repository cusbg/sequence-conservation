FROM Debian:buster-slim

RUN apt-get update

RUN apt-get install wget python2 python3 \
 && apt-get install make g++ zlib1g-dev

WORKDIR /opt/conservation

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz \
 && wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz \
 && wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz \
 && wget http://compbio.cs.princeton.edu/conservation/conservation_code.tar.gz \
 && tar -xf ncbi-blast-2.9.0+-x64-linux.tar.gz \
 && tar -xf cd-hit-v4.8.1-2019-0228.tar.gz \
 && tar -xf muscle3.8.31_i86linux64.tar.gz \
 && tar -xf conservation_code.tar.gz \
 && cd cd-hit-v4.8.1-2019-0228 \
 && make \
 && cd ..
 && rm *.gz

ENV BLASTDB /data/conservation/blast-database

VOLUME ["/data/conservation"]

COPY calculate_conservation.py ./
COPY prepare_database.sh ./

CMD ["/bin/bash"]
