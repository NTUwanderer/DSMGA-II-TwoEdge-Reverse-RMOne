#!/bin/bash
if [ "$#" == "0" ]; then
    echo "./sw-all.sh DIR"
    exit 1
fi
DIR=$1
mkdir -p $DIR

nohup nice -n19 ./sweep 100 10 1 > $DIR/mktrap_100 &
nohup nice -n19 ./sweep 120 10 2 > $DIR/folded_120 &
nohup nice -n19 ./sweep 100 10 3 > $DIR/cyclic_100 &
nohup nice -n19 ./sweep-sat.sh $DIR/SAT 100
nohup nice -n19 ./sweep-spin.sh $DIR/SPIN 400


