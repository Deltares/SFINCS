#!/bin/sh
# add autoconf
# Bug https://bugzilla.redhat.com/show_bug.cgi?id=1869030
# export CONFIG_SHELL=/bin/bash
./autogen.sh
# configure software
FCFLAGS=-fallow-argument-mismatch ./configure
# build
make
# install
make install
