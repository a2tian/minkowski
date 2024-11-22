#!/bin/bash

# Install dependencies
yum install -y gmp gmp-devel mpfr mpfr-devel autoconf automake libtool pkgconfig

# Install fplll
wget https://github.com/fplll/fplll/releases/download/5.4.5/fplll-5.4.5.tar.gz
gunzip fplll-5.4.5.tar.gz
tar xf fplll-5.4.5.tar
cd fplll-5.4.5
./autogen.sh
./configure --without-qd
make
make install
cd ..
rm -rf fplll-5.4.5

# Install NTL
wget https://libntl.org/ntl-11.5.1.tar.gz
gunzip ntl-11.5.1.tar.gz
tar xf ntl-11.5.1.tar
cd ntl-11.5.1/src
./configure SHARED=on
make
make install
cd ../..
rm -rf ntl-11.5.1
