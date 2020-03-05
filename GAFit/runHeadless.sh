#!/bin/bash

touch ./headlessOut/$1.lck
echo $1 > ./headlessOut/$1.log
matlab -nosplash -nodesktop -nojvm -r "run('main.m');quit ">./headlessOut/$1.log
rm ./headlessOut/$1.lck

