#!/bin/bash -x

sudo add-apt-repository ppa:cae-team
sudo apt-get update

# install development packages
sudo apt-get install dpkg-dev git core

#
# install dependecies for OpenFOAM
#
sudo apt-get install cdbs debhelper binutils-dev flex  libparmetis-dev libreadline-dev zlib1g-dev doxygen graphviz texlive-base  openmpi-bin libscotch-dev libxt-dev

sudo apt-get install libmgridgen-dev libmesquite-dev

sudo apt-get upgrade
