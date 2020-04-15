#!/bin/bash

IN_DIR=${1}
OUTPUT_DIR=${2}
CWD=${2}

declare -x KB_FLAG=absent

mkdir -p ${OUTPUT_DIR}

declare -x OUTPUT_FN="${OUTPUT_DIR}/summary.txt"

rm -f ${OUTPUT_FN}
echo -e "chr\tk\tflag\trealtime\tusertime\tsystime" >> ${OUTPUT_FN}

for K in `seq 2 12`
do
  for SUFFIX in `seq 1 22` X Y
  do
    CHROMOSOME="chr${SUFFIX}"
    echo "${CHROMOSOME} | ${K} | ${KB_FLAG}"
    declare -x IN_FN="${IN_DIR}/${CHROMOSOME}.${K}.${KB_FLAG}.stderr.log"
    echo "${IN_FN}"
    awk -v OFS="\t" -v CHR=${CHROMOSOME} -v K=${K} -v KB_FLAG=${KB_FLAG} '{ print CHR, K, KB_FLAG, $2, $4, $6 }' ${IN_FN} >> ${OUTPUT_FN}
  done
done