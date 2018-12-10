#!/bin/bash
####################################
# User input configuration 
# intput: 	$1: scannerID
#		$2: config.txt
#		duration: X in seconds
#		imagepath: example /home/spaldinglab/Desktop/images
#		scannerID: (RootHub in 3 digits)_(port number in a single digit)
#		scannerID example: root hub-2, port-4 ==> 002_4
####################################
####################################
##############################################
#machType="Laptop"
filePath="/home/spalding/Desktop/ScannerAuto"
#changes 'date' format by replacing ' ' to '_' and ':' to '-'
#1)It is not good to have ' ' in the file/folder name
#2)snapper does not accept ':' in the file/folder name
status_name=$(echo "$(date)" | sed -e 's/ /_/g'| sed -e 's/:/-/g')
#name for status text file
status_name="status_file_"$status_name
#name for log file
stdout_name="stdout_file_"$status_name
#takes logs into text file? and call 'run_scanner.sh'
$filePath/run_scanner.sh $1 $status_name < $2 1 > $stdout_name.txt $3
