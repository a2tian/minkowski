#!/bin/bash

# Install dependencies
yum install -y gmp gmp-devel mpfr mpfr-devel autoconf automake 

# Install fplll
git clone https://github.com/fplll/fplll.git
cd fplll
./autogen.sh
./configure
make
make install
cd ..

# Install NTL
wget https://libntl.org/ntl-11.5.1.tar.gz
gunzip ntl-11.5.1.tar.gz
tar xf ntl-11.5.1.tar
cd ntl-11.5.1/src
./configure
make
make install
cd ../..