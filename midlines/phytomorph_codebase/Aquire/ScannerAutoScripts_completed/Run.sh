#!/bin/bash
# Run scanner_api.sh only with scanner Id input
# Find computer number and store it then pass it to scanner_api.sh to look up ScannerWiring.txt
## April.8.2015.
### April.14.2015.
####################################
# User input configuration 
# intput: 	$1: scannerID 
####################################
####################################

## Get IP address fot this machine
IP=$(ifconfig | grep 'inet addr' | head -1 | sed -e 's/ /+/g' | cut -d ':' -f2 |cut -d '+' -f1)
echo 'Computer IP address requested is: ' $IP
## Get its computer number from the ip table
this_computer=$(cat /mnt/snapper/kernelSwellingData/config/IPTable.txt | grep $IP | cut -d ':' -f1)
echo 'The scanner is running on computer: ' $this_computer

## store ScannerWiring.txt into line
#line=$(cat $path/ScannerWiring.txt)
## pass $this_computer to scanner_api.sh to look for device number to pass to run_scanner.sh
#second half of the line after '$1', is hard coded 'config.txt' input. '&' allows run this as a background.
/home/spaldings_lab/Desktop/Scripts/Scripts_main/scanner_api.sh $1 /home/spaldings_lab/Desktop/Scripts/Scripts_main/new_config.txt $this_computer &



## I need to redirect stdout file
## Confusion in passing variables
## $this_computer does not seem to be pssed


####I can scan at 400 or 800 or  resolution but somehow I cannot scan at 100, 1200, 1600, 2400, 3200, 4800.
####Error message comes up
