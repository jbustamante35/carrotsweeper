#!/bin/bash
# kill all scanner
killall -9 Scott_run.sh
killall -9 run_single_scan
killall -9 scanner_api.sh
killall -9 run_scanner.sh
killall -9 'sleep'

#pkill -9 <process_name> It does not work in shell script
#killall -9 <process_name> It does work in Shell script


#pgrep <expression> This extract only number PID








# store process ID for running scanners
#In which case do I get a space?
#all_scanner=$(ps | grep run_scanner.sh | cut -d" " -f1)
# I need to count number of digits of process ID
# add space at the end to ease the count
#echo $$all_scanner | wc -c

#all_scanner=$all_scanner" "

#echo $(($(echo $$all_scanner | wc -c)))
# count total number of running scanners
#all_scanner_count=$(($(echo $$all_scanner | wc -c)/2))

#echo $all_scanner_count
# if counter is not 0
#if [ $all_scanner_count != 0 ]; then
# for i is 1:counter
#for i in $all_scanner_count; do
# kill all the running scanners
#kill -9 $all_scanner | cut -d" " -f$i
#((all_scanner_count++))
#done
#fi

#echo $all_scanner_count

# if no scanner is running
#if [ $all_scanner_count == 9 ]; then
#echo "Currently no scanners to kill"
#fi

#deal with when no scanner is running




















# store currently scanning scanners as epson2:libusb:(root hub):(device) 
#current=$(ps ux | grep "epson2:libusb" | grep -v "/time" | cut -d "-" -f4 | cut -d " " -f2)
# display currently scanning scanners
#echo "These are currently running scanner is: " $current
# counting number of scanners those are slpeeing and store
#sleepingcount=$(ps ux | grep "sleep" | grep -v "grep" | grep "sleep" -c)
# display number of scanners those are sleeping
## I may add funtion that can change the string is/are corresponding to the situation
#echo "The number of sleeping scanners are: " $sleepingcount
#sleepingtime=$(ps ux | grep "sleep" | grep -v "grep" | grep "sleep")
# chop the variable that contains long string/ or possibly multiple strings
#newsleepingtime=$(ps ux | grep "sleep" | grep -v "grep" | grep "sleep" | cut -d "S" -f2 | cut -d ":" -f3 | cut -d "0" -f3)
#echo "Sleeping time left are: "$newsleepingtime

############see process and kill process
#Ubuntu has a command that lists all current processes.If you find the process ID,you can kill it by kill command.
#Type following to list all your own running processes
#ps ux
#And use command to kill any process you want,the second column is PIDnumber.
#kill -9 PIDnumber
#Ubuntu also has a graphical program to view current processes.Navigate to System–>Administrator–>System Monitor,you will see system processes #at Process tab,kill any by select it and click “End Process” button.
#Name my process with something I know or want 



##############################################
# 
##############################################
## Feb18, how can I attach scannerID into log?
## maybe we need to modify run_scanners.sh or run_single_scanner.sh to leave mark for scanner-identification on folder..?
## I only can show the status while the scanner is actually running. I cannot track if the scanner reveived order or not while 
## the scanner is sleeping.
## Add renaming?
## use grep and cut to show multiple strings stored in the variable nicely and avoide problem with strings with multiple stars
