#!/bin/bash
####################################
# this is general downloader for ScannerAuto
# default location for program is /home/
# vars intput:
#	$1: fileloc: file location
#	$2: repo: repository address
#	$3: filename
####################################



# download [filename]

echo "start download: $3 from $2"
sudo wget -O $1/$3 $2/$3
sudo chmod +x $1/$3
echo "end download: $3 from $2"

