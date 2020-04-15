#!/bin/bash

KB_BIN=${1}
K=${2}
FLAG=${3}
IN_FN=${4}
OUT_FN=${5}

TIMEFORMAT='real: %3R user: %3U sys: %3S cpu: %P'
time ${KB_BIN} --k=${K} --${FLAG} < ${IN_FN} > ${OUT_FN}