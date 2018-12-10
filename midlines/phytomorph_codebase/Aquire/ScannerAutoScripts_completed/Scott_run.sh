#!/bin/bash
#command to run scanner_api.sh only with scanner Id input
####################################
# User input configuration 
# intput: 	$1: scannerID
####################################
####################################
#second half of the line after '$1', is hard coded 'config.txt' input. '&' allows run this as a background.
/home/spaldinglab/Desktop/production_scripts/Scripts_main/scanner_api.sh $1 /home/spaldinglab/Desktop/production_scripts/Scripts_main/config.txt &
