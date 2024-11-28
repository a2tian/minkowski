#!/bin/bash

# Install dependencies
yum install -y gmp gmp-devel mpfr mpfr-devel autoconf automake libtool pkgconfig

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
