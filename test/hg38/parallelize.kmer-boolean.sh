#!/bin/bash

KB_BIN=${1}
OUTPUT_DIR=${2}
FASTA_DIR=${3}
CWD=${4}

declare -x SLURM_PARTITION="queue0"
declare -x SLURM_MEM_PER_CPU="1G"
declare -x SLURM_DESCRIPTION="kb"
declare -x SLURM_PROCESSORS=1

declare -x KB_FLAG=absent

mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/logs

for K in `seq 2 12`
do
  for SUFFIX in `seq 1 22` X Y
  do
    CHROMOSOME="chr${SUFFIX}"
    echo "${CHROMOSOME} | ${K} | ${KB_FLAG}"
    declare -x IN_FN="${FASTA_DIR}/${CHROMOSOME}.fa"
    declare -x OUT_FN="${OUTPUT_DIR}/${CHROMOSOME}.${K}.${KB_FLAG}.txt"
    declare -x SLURM_CMD="${CWD}/slurm.kmer-boolean.sh ${KB_BIN} ${K} ${KB_FLAG} ${IN_FN} ${OUT_FN}"
    sbatch \
      --parsable \
      --partition=${SLURM_PARTITION} \
      --job-name="${SLURM_DESCRIPTION}-${CHROMOSOME}-${K}" \
      --nodes=1 \
      --ntasks-per-node=${SLURM_PROCESSORS} \
      --mem-per-cpu=${SLURM_MEM_PER_CPU} \
      --output="${OUTPUT_DIR}/logs/${CHROMOSOME}.${K}.${KB_FLAG}.stdout.log" \
      --error="${OUTPUT_DIR}/logs/${CHROMOSOME}.${K}.${KB_FLAG}.stderr.log" \
      --wrap="${SLURM_CMD}"
  done
done
