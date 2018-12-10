#!/bin/bash
#command to run scanner_api.sh only with scanner Id input
####################################
# User input configuration 
# intput: 	$1: scannerID
####################################
####################################
filePath="/home/spalding/Desktop/ScannerAuto"
#second half of the line after '$1', is hard coded 'config.txt' input. '&' allows run this as a background.
$filePath/scanner_api.sh $1 /$filePath/config.txt $2 & 
