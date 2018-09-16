FROM informaticsmatters/rdkit_debian_base
MAINTAINER Tim Dudgeon <tdudgeon@informaticsmatters.com>

RUN apt-get update &&\
 apt-get upgrade -y &&\
 apt-get clean -y

RUN apt-get install git libeigen3-dev gfortran unzip libxml2-dev zlib1g-dev libblas-dev liblapack-dev vim -y
RUN apt-get install cmake locate -y
RUN apt-get install libf2c2-dev -y

RUN updatedb
RUN pwd
RUN ls
RUN locate libblas
RUN git clone https://github.com/thegodone/yaehmop.git
WORKDIR yaehmop/tightbind

RUN make

RUN cd /yaehmop/tightbind
RUN ls
RUN ./bind example.bind eht_parms.dat
