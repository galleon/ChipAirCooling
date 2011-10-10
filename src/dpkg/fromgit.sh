#!/bin/bash -x

rm -rf OpenFOAM-1.6-ext openfoam-1.6-ext_*
git clone git://openfoam-extend.git.sourceforge.net/gitroot/openfoam-extend/OpenFOAM-1.6-ext
#git clone Master_OpenFOAM-1.6-ext OpenFOAM-1.6-ext
cd OpenFOAM-1.6-ext

git branch -D eads
git checkout -b eads 8673d8a
git merge origin/packaging/ubuntu/10.04

mkdir -p tutorials/heatTransfer/chipAirCoolingFoam
cp -a ../../../run/* tutorials/heatTransfer/chipAirCoolingFoam 
cp -a ../../*Foam applications/solvers/heatTransfer
cp -a ../../airbusMaterialModels src

patch -p1 < ../patchAllwmake

git add -A
git commit -a -m "Copied chipAirCooling project"

dpkg-buildpackage 2>&1 | tee ~/dpkg-buildpackage.log
