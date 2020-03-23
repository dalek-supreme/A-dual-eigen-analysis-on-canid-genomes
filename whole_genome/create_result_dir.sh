#!/bin/bash
cd ~
mkdir enrichment
cd enrichment
for layer in {1..$NUM_LAYERS}
    do
        mkdir $layer
        cd $layer
        mkdir positive
        mkdir negative
        cd -
    done
