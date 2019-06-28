#!/bin/sh

seq 1 4 | xargs -n 1 -P 4 bash Maud_batch.sh

