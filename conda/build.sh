#!/bin/sh

cp -r $SRC_DIR/src/*.py $PREFIX/bin
#mkdir $PREFIX/models
#cp -r $SRC_DIR/models/* $PREFIX/models

#Create a symlink between ffpe_filtering and cli.py (command line interface) script 
ln -s $PREFIX/bin/cli.py $PREFIX/bin/ffpe_filtering
#This makes the symlink executable
chmod +x $PREFIX/bin/ffpe_filtering
