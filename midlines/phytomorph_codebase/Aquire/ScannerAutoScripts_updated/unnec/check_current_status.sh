#!/bin/bash
####################################
# given device number( or root hub and port), 
# check which scanners are currently running, what image number the scanner is currently working on and such. 
####################################
# User input configuration 
####################################
####################################
# grep the scanner ID that is associated with date above 
numscanners=$(ps aux | grep run_scanner | grep -v color=auto | cut -d "." -f4 | cut -d " " -f2| sed -e 's/ /_/g')
numscanners=$numscanners" "

echo $numscanners

# count number of scanners by counting the number of spaces
counter=$(($(echo $numscanners | wc -c)/2))
echo $counter 
if [ $counter != 0 ]; then 
#while [ $counter != 0 ]; do
i=1
for i in $counter; do
current_scan_is[i]=$(echo $numscanners | cut -d" " -f$i)
#current_scannerID=$(ps aux | grep run_scanner | grep -v color=auto | cut -d "." -f4 | cut -d " " -f2| sed -e 's/ /_/g')
scanner_status_file[i]=$(ps aux | grep run_scanner | grep -v color=auto | cut -d "." -f4 | cut -d " " -f3 | cut -d " " -f3  | cut -d " " -f$i)
#sed 's\:.*//'
# take .log lines into variable
scanner_status_file[i]=$1${scanner_status_file[i]}".log" 
echo ${scanner_status_file[i]}
# contains only the last line of scanner_status_file
scanner_status_show[i]=$(cat ${scanner_status_file[i]} | tail -1)

# divive scanner_status_show into two pieces
show_part1[i]=$(echo ${scanner_status_show[i]} | cut -d: -f1)
show_part2[i]=$(echo ${scanner_status_show[i]} | cut -d: -f2)
# this displays on the screen
echo ${show_part1[i]}":" "Scanner number "${current_scan_is[i]}" is" ${show_part2[i]}"."
((i++))
done
fi

if [ $counter == 0 ]; then
echo "No scanner is running"
fi


