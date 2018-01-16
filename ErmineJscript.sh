#!/usr/bin/bash

set -u
#cd "$HOME"/"$ERMINEJ_HOME"

scoreFile=$1
scoreCol=$2
annoFile=$3
outFile=$4

./ermineJ.sh --scoreFile $scoreFile --classFile ~/ermineJ.data/go_daily-termdb.rdf-xml.gz --annots $annoFile --test GSR --stats PRECISIONRECALL  --scoreCol $scoreCol -b  --minClassSize 10 > $outFile
