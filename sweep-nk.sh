#!/bin/bash
REPEAT=99
nThreads=$((`nproc --all` - 0))
if [ "$#" == "0" ]; then
echo "./sweep-sat step DIR ell1 ell2 ... ellN"
exit 1
fi
step=$1
DIR=$2
mkdir -p $DIR
while [ $# -gt 2 ]
do
    for (( i=0; $i<=$REPEAT; i=$i+1 ))
    do
        NUM=`echo $i | awk '{printf "%03d", $1}'`
        ELL=$3
        nohup nice -n19 ./sweep $ELL 10 4 $step $i > ./$DIR/$ELL-$NUM &
        echo "Submitting $ELL-$NUM"
        sleep 2
        TT=$(ps -xaf | grep sweep | wc -l)
        while [ $TT -gt $nThreads ]
        do
            sleep 1
            TT=$(ps -xaf | grep sweep | wc -l)
        done
    done
    shift
done
