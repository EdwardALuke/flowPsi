#!/bin/bash
set -e
set -u

INSTALL_DIR=${FLOWPSI_INSTALL_PATH-/usr/local}

INSTALL_PATH=$INSTALL_DIR/$FLOWPSI_INSTALL_DIR

echo INSTALL_PATH = $INSTALL_PATH

echo Making Directories
mkdir -p $INSTALL_PATH
mkdir -p $INSTALL_PATH/lib
mkdir -p $INSTALL_PATH/bin
mkdir -p $INSTALL_PATH/docs

cp flowpsi.conf $INSTALL_PATH/flowpsi.conf
cp revision.conf $INSTALL_PATH/revision.conf

echo Installing Library Files
LIB_POSTFIX="so"

ARCH=${LOCI_ARCH-`uname -s`}
if [ $ARCH == "Darwin" ]; then
    LIB_POSTFIX="dylib"
fi

echo Installing flowpsi
cp src/flowpsi $INSTALL_PATH/bin

cp lib/* $INSTALL_PATH/lib

mkdir -p $INSTALL_PATH/docs
if [ -e guide/guide.pdf ] ; then
  cp guide/guide.pdf $INSTALL_PATH/docs/flowPsiGuide.pdf
fi


echo Installing \#include files
mkdir -p $INSTALL_PATH/include
cp include/*.h include/*.lh $INSTALL_PATH/include

# make directories open to everyone
chmod -R a+rX $INSTALL_PATH
