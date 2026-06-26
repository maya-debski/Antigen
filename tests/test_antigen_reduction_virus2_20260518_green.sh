#!/usr/bin/env bash
#---------------------------------------------------------
# Test data reduction, Antigen, VIRUS2, 2026-05-18 data files, just the GREEN channel
# UT_BOX_DATA_PATH
#     "UT BOX / Greg Zeimann Share / data / VIRUS2 / 20260518"
#---------------------------------------------------------

ANTIGEN_ROOT=$HOME/Code/github-com/mcdo-hjst/Antigen
DATA_DATE=20260518
DATA_ROOT=${ANTIGEN_ROOT}/wip/test_data/test_reduction_virus2_data_${DATA_DATE}
REDUCED_ROOT=${DATA_ROOT}/reduced/VIRUS2/${DATA_DATE}

echo "ANTIGEN_ROOT = ${ANTIGEN_ROOT}"
echo "DATA_DATE    = ${DATA_DATE}"
echo "DATA_ROOT    = ${DATA_ROOT}"
echo "REDUCED_ROOT = ${REDUCED_ROOT}"

antigen_base_reduction.py -i ${DATA_ROOT} -o ${REDUCED_ROOT} -c ${DATA_DATE} -m VIRUS2 -v -d -u D2G

antigen_make_cubes.py -r ${REDUCED_ROOT}/VIRUS2_D2G -o ${REDUCED_ROOT}/VIRUS2_D2G/cubes -v -d
