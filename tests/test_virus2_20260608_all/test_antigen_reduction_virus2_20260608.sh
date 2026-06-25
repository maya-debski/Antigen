#!/usr/bin/env bash
#---------------------------------------------------------
# Test data reduction, Antigen, VIRUS2, 2026-06-08 data files
# UT_BOX_DATA_PATH
#     "UT BOX / VIRUS2 Folder / HJST_SHARE / 2026_engineering_run_datasets / oppor_20260608_night /"
#---------------------------------------------------------

ANTIGEN_ROOT=$HOME/Code/github-com/mcdo-hjst/Antigen
TEST_DATA_DATE=20260608
TEST_DATA_ROOT=${ANTIGEN_ROOT}/wip/test_data/test_reduction_virus2_data_20260608/oppor_20260608_night
TEST_REDUCED_ROOT=${TEST_DATA_ROOT}/reduced/VIRUS2/${TEST_DATA_DATE}

echo "ANTIGEN_ROOT      = ${ANTIGEN_ROOT}"
echo "TEST_DATA_DATE    = ${TEST_DATA_DATE}"
echo "TEST_DATA_ROOT    = ${TEST_DATA_ROOT}"
echo "TEST_REDUCED_ROOT = ${TEST_REDUCED_ROOT}"

CHANNEL_IDS=("D2D" "D2R" "D2G" "D2B")
for CHAN_ID in "${CHANNEL_IDS[@]}"; do
    echo "----------------------------------------"
    echo "Processing channel: ${CHAN_ID}, base reduction"
    antigen_base_reduction.py -i ${TEST_DATA_ROOT} -o ${TEST_REDUCED_ROOT} -c ${TEST_DATA_DATE} -m VIRUS2 -v -d -u ${CHAN_ID} -w 7

    echo "----------------------------------------"
    echo "Processing channel: ${CHAN_ID}, make cube"
    antigen_make_cubes.py -r ${TEST_REDUCED_ROOT}/VIRUS2_${CHAN_ID} -o ${TEST_REDUCED_ROOT}/VIRUS2_${CHAN_ID}/cubes -v -d
done

echo "Channel processing complete"
