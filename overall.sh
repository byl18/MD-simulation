#!/bin/bash
# trails=("trail02" "trail03" "trail04" "trail05" "trail06" "trail07" "trail08" "trail09" "trail10")
trails=("trail02")
for trail in "${trails[@]}"; do

echo "Entering: $trail"
cd "$trail" || { echo "Dir not exist $trail"; exit 1; }
echo "In $trail running md_total.sh..."
nohup bash md_total.sh &
PID=$!
echo "process PID: $PID"
wait $PID
echo "$trail done"
cd ..
echo "----------------------------------------"

done
echo "All trails are done！"


