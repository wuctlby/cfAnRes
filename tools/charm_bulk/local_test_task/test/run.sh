OPTION="-b --configuration json://cfg.json"
MEMORY="--aod-memory-rate-limit 4000000000 --shm-segment-size 24000000000"
# CPU="--pipeline hf-correlator-charm-hadrons-reduced:1"
# export OMP_NUM_THREADS=8
LOGFILE="log.txt"

input_path=$(head -n 1 input_data.txt)
dir=$(dirname "$input_path")
# parent_dir=$(realpath "$dir/parent")
local_dir="/home/wuct/MetaData/DATA/PbPb/2023/pass4/Meta/alice"
parent_dir="/home/wuct/MetaData/DATA/PbPb/2023/pass4/Meta/parent"


echo "Input path: $input_path"
echo "Parent directory: $parent_dir"

o2-analysis-hf-pid-creator $OPTION |\
o2-analysis-ft0-corrected-table $OPTION |\
o2-analysis-pid-tof-base $OPTION |\
o2-analysis-pid-tof-full $OPTION |\
o2-analysis-tracks-extra-v002-converter $OPTION |\
o2-analysis-multcenttable $OPTION |\
o2-analysis-event-selection-service $OPTION |\
o2-analysis-propagationservice $OPTION |\
o2-analysis-pid-tpc-service $OPTION |\
o2-analysis-trackselection $OPTION |\
o2-analysis-hf-candidate-selector-d0 $OPTION |\
o2-analysis-hf-candidate-creator-2prong $OPTION |\
o2-analysis-hf-task-pt-fluc-charm-hadrons $OPTION --aod-file @input_data.txt --aod-parent-access-level 1 --aod-parent-base-path-replacement "alien://;$parent_dir" $MEMORY > $LOGFILE 2>&1