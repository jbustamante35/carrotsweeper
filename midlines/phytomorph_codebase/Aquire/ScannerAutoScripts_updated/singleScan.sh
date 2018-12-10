#!/bin/bash
####################################
# run single scanner and save single image
# vars intput: 	$1: file name
#		$2: resolution
#		$3: root hub
#		$4: device number
#		$5: x 
#		$6: y
####################################

echo "Enter x in mm from top-right"
read x_crop
echo "Enter y in mm from top-right"
read y_crop
echo "Enter resolution"
read res
echo "Enter directory"
read tory
echo "Enter filename with no space"
read filename
folder="/home/spaldinglab/Documents/$tory"
mkdir $folder
echo "Your image is in $folder" 
scannerIDdevnum=$(/home/spaldinglab/Desktop/production_scripts/Scripts_main/getDeviceForPort.sh 002 3)
scanTime=$(/usr/bin/time -f %e 2>&1 scanimage --mode color -x $x_crop -y $y_crop -d "epson2:libusb:002:$scannerIDdevnum" --resolution $res --format tiff > $folder/$filename.tif)

