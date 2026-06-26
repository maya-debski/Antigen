#!/usr/bin/env bash
#---------------------------------------------------------
# Test data reduction
# Antigen, VIRUS2, 2026-06-08, GREEN channel only
# Data source: "UT BOX / VIRUS2 Folder
#                      / HJST_SHARE
#                      / 2026_engineering_run_datasets
#                      / oppor_20260608_night "
#---------------------------------------------------------

ANTIGEN_ROOT=$HOME/Code/github-com/mcdo-hjst/Antigen
DATA_DATE=20260518
DATA_ROOT=${ANTIGEN_ROOT}/tests/oppor_20260608_night
REDUCED_ROOT=${DATA_ROOT}/reduced/VIRUS2/${DATA_DATE}
MANIFEST_FILE=${ANTIGEN_ROOT}/tests/manifest_base_reduction_minimal_vestuto.yml
CHAN_ID=D2G

echo "----------------------------------------"
echo "Bash script variables:"
echo " CHAN_ID      = ${CHAN_ID}"
echo " DATA_DATE    = ${DATA_DATE}"
echo " ANTIGEN_ROOT = ${ANTIGEN_ROOT}"
echo " DATA_ROOT    = ${DATA_ROOT}"
echo " REDUCED_ROOT = ${REDUCED_ROOT}"

echo "----------------------------------------"
echo "Running antigen_base_reduction.py:"
antigen_base_reduction.py \
    -i ${DATA_ROOT} \
    -o ${REDUCED_ROOT} \
    -c ${DATA_DATE} \
    -m VIRUS2 \
    -M ${MANIFEST_FILE} \
    -u ${CHAN_ID} \
    -w 7 \
    -j 1 \
    -v \
    -d \

echo "----------------------------------------"
echo "Running antigen_make_cubes.py:"
antigen_make_cubes.py \
    -r ${REDUCED_ROOT}/VIRUS2_${CHAN_ID} \
    -o ${REDUCED_ROOT}/VIRUS2_${CHAN_ID}/cubes \
    -v \
    -d \

echo "----------------------------------------"
echo "Done"
