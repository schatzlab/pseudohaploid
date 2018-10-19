#!/bin/sh

## Minimum alignment Identity
MIN_IDENTITY=.90

## Minimum alignment length to consider
MIN_LENGTH=1000

GENOME=$1

nucmer --maxmatch -c 100 -l 500 $GENOME $GENOME -p self
delta-filter -l $MIN_LENGTH -i $MIN_IDENTITY self.delta > self.filter.delta

$BIN/pseudohaploid.chain.pl self.filter.delta > self.contained
grep '^#' self.contained > self.contained.id
$BIN/filter_seq -v self.contained.ids $GENOME > filtered.fa
