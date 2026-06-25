#!/usr/bin/env bash
#---------------------------------------------------------
# Test data reduction
# Antigen, VIRUS2, 2026-06-08, DRED channel only
# Data source: "UT BOX / VIRUS2 Folder / HJST_SHARE
#                 / 2026_engineering_run_datasets / oppor_20260608_night "
#---------------------------------------------------------

ANTIGEN_ROOT=$HOME/Code/github-com/mcdo-hjst/Antigen
DATA_DATE=20260608
DATA_ROOT=${ANTIGEN_ROOT}/tests/test_data_virus2_20260608
REDUCED_ROOT=./reduced
MANIFEST_FILE=./manifest.yml
CHAN_ID=D2D

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
    --obs_date ${DATA_DATE} \
    --instrument VIRUS2 \
    --manifest ${MANIFEST_FILE} \
    --unit-id ${CHAN_ID} \
    --time_radius 7 \
    --good_arc_residual_limit 1 \
    --recipe base_reduction \
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
