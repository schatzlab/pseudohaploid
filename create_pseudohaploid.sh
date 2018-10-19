#!/bin/sh

set -o pipefail
set -xv

BIN="$( cd "$(dirname "$0")" ; pwd -P )"

## Minimum alignment Identity
MIN_IDENTITY=90

## Minimum alignment length to consider
MIN_LENGTH=1000

## Minimum containment percentage from overlap chains to filter contig
MIN_CONTAIN=93

if [ $# -lt 2 ]
then
  echo "create_pseudohaploid.sh assembly.fa outprefix"
  exit
fi

GENOME=$1
PREFIX=$2

## You may want to replace this with sge_mummer for large genomes
## See: https://github.com/fritzsedlazeck/sge_mummer
nucmer --maxmatch -c 100 -l 500 $GENOME $GENOME -p $PREFIX

## Pre-filter for just longer high identity alignments
delta-filter -l $MIN_LENGTH -i $MIN_IDENTITY $PREFIX.delta > $PREFIX.filter.delta

## Create the coord file
show-coords -rclH $PREFIX.filter.delta > $PREFIX.filter.coords

## Find and analyze the alignment chains
## Note you can rerun this step multiple times from the same coords file
$BIN/pseudohaploid.chains.pl $PREFIX.filter.coords $MIN_IDENTITY $MIN_CONTAIN > $PREFIX.chains

## Generate a list of contained contigs
## This can also be rerun from the same chain file but using different containment thresholds
grep '^#' $PREFIX.chains | \
  awk -v cut=$MIN_CONTAIN '{if ($4 >= cut){print ">"$2}}' > $PREFIX.contained.ids

## Finally filter the original assembly to remove the contained contigs
$BIN/filter_seq -v $PREFIX.contained.ids $GENOME > $PREFIX.pseudohap.fa
