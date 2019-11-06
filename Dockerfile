FROM debian:buster-slim

ARG USER=5988

RUN apt-get update

RUN apt-get -y install wget curl python2 python3 \
 && apt-get -y install make g++ zlib1g-dev \
 && apt-get -y install python-numpy

RUN addgroup --gid ${USER} user \
 && useradd --uid ${USER} --gid ${USER} user

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
 && cd .. \
 && rm *.gz

ENV BLASTDB /data/conservation/blast-database

COPY prepare_database.sh ./
COPY download_database.sh ./
COPY calculate_conservation.py ./

RUN chown user:user -R /opt/conservation \
 && chmod a+x calculate_conservation.py \
 && chmod a+x download_database.sh \
 && chmod a+x prepare_database.sh \
 && mkdir -p /data/conservation \
 && chown user:user /data/conservation

USER ${USER}:${USER}

VOLUME ["/data/conservation"]

CMD ["/bin/bash"]
