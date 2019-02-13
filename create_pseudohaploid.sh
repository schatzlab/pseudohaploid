#!/bin/sh

set -o pipefail

BIN="$( cd "$(dirname "$0")" ; pwd -P )"

## Minimum alignment Identity
MIN_IDENTITY=90

## Minimum alignment length to consider
MIN_LENGTH=1000

## Minimum containment percentage from overlap chains to filter contig
MIN_CONTAIN=93

## Maximum distance in bp allowed between alignments on the same alignment chain
MAX_CHAIN_GAP=20000

if [ $# -lt 2 ]
then
  echo "create_pseudohaploid.sh assembly.fa outprefix"
  exit
fi

GENOME=$1
PREFIX=$2

echo "Generating pseudohaploid genome sequence"
echo "----------------------------------------"
echo "GENOME: $GENOME"
echo "OUTPREFIX: $PREFIX"
echo "MIN_IDENTITY: $MIN_IDENTITY"
echo "MIN_LENGTH: $MIN_LENGTH"
echo "MIN_CONTAIN: $MIN_CONTAIN"
echo "MAX_CHAIN_GAP: $MAX_CHAIN_GAP"
echo

## You may want to replace this with sge_mummer for large genomes
## See: https://github.com/fritzsedlazeck/sge_mummer
echo "1. Aligning $GENOME to itself with nucmer"
(nucmer --maxmatch -c 100 -l 500 $GENOME $GENOME -p $PREFIX) >& nucmer.log
numorig=`grep -c '^>' $GENOME`
echo "Original assembly has $numorig contigs"
echo

## Pre-filter for just longer high identity alignments
echo "2. Filter for alignments longer than $MIN_LENGTH bp and below $MIN_IDENTITY identity"
delta-filter -l $MIN_LENGTH -i $MIN_IDENTITY $PREFIX.delta > $PREFIX.filter.delta
echo

## Create the coord file
echo "3. Generating coords file"
show-coords -rclH $PREFIX.filter.delta > $PREFIX.filter.coords
echo

## Find and analyze the alignment chains
## Note you can rerun this step multiple times from the same coords file
echo "4. Identifying alignment chains: min_id: $MIN_IDENTITY min_contain: $MIN_CONTAIN max_gap: $MAX_CHAIN_GAP"
($BIN/pseudohaploid.chains.pl $PREFIX.filter.coords \
  $MIN_IDENTITY $MIN_CONTAIN $MAX_CHAIN_GAP > $PREFIX.chains) >& $PREFIX.chains.log
cat $PREFIX.chains.log
echo

## Generate a list of contained contigs
## This can also be rerun from the same chain file but using different containment thresholds
echo "5. Generating a list of redundant contig ids using min_contain: $MIN_CONTAIN"
grep '^#' $PREFIX.chains | \
  awk -v cut=$MIN_CONTAIN '{if ($4 >= cut){print ">"$2}}' > $PREFIX.contained.ids
numcontained=`wc -l $PREFIX.contained.ids | awk '{print $1}'`
echo "Identified $numcontained redundant contig to remove in $PREFIX.contained.ids"
echo


## Finally filter the original assembly to remove the contained contigs
echo "6. Creating final pseudohaploid assembly in $PREFIX.pseudohap.fa"
($BIN/filter_seq -v $PREFIX.contained.ids $GENOME > $PREFIX.pseudohap.fa) >& $PREFIX.pseudohap.fa.log
numfinal=`grep -c '^>' $PREFIX.pseudohap.fa`
echo "Pseudohaploid assembly has $numfinal contigs"
