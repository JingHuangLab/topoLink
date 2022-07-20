#!/bin/bash

### You need to create your conda environment
# conda activate topolinks


## parameters
inPath=$1  ## Path contain files to analysis
outPath=${2:-"/tmp/topo_links"} ## Path to save the files of results
summaryfile=${3:-'summary_topo_links.txt'}
detail=${4:-'no-detail'} ## or 'detail'
scanend=${5:-36}
scanbegin=${6:-4}
threshscore=${7:-0.8}
discutoff=${8:-10}
ter_rm=${9:-15}

filelist=$(ls $inPath | grep -w pdb)' '$(ls $inPath | grep -w cif)

for file in $filelist
do
    infile=$inPath'/'$file
    echo $file
    python topo_links.py  -in $infile  -out $outPath \
        -sf $summaryfile  --$detail  -se $scanend  -sb $scanbegin  \
        -ts $threshscore  -d $discutoff -rm $ter_rm
done
