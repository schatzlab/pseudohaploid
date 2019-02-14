#!/bin/sh

echo "Running the simple example"
echo "==============================================================="
mkdir -p simple
cd simple
ln -sf ../simple.fa
../../create_pseudohaploid.sh simple.fa ph.simple
cd ..

echo
echo "This should report: Pseudohaploid assembly has 1 contigs"

echo "==============================================================="

echo
echo


echo "Running the basic example"
echo "==============================================================="
mkdir -p basic
cd basic
ln -sf ../basic.fa
../../create_pseudohaploid.sh basic.fa ph.basic
cd ..

echo
echo "This should report: Pseudohaploid assembly has 3 contigs"
echo "==============================================================="
