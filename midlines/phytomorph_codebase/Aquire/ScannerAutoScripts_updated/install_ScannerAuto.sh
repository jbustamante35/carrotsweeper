#!/bin/bash
####################################
# this is general installation for ScannerAuto
# default location for program is /home/
# vars intput:
#	
####################################
repo="https://bitbucket.org/[username]/[reponame]/raw"
fileloc="/home/spalding/Desktop/ScannerAuto"
sudo mkdir -p $fileloc/

echo "Do you want to Download+Install ScannerAuto, Install ScannerAuto or EXIT?"
select selection in "DownNInstall" "InstallOnly" "EXIT"; do
case $selection in
	DownNInstall )
		sudo wget -O $fileloc/downloader.sh $repo/downloader.sh;
		sudo chmod +x $fileloc/downloader.sh;
		sudo $fileloc/downloader.sh $fileloc $repo [filename];
		break;;
	InstallOnly )
		break;;
	EXIT )
		exit;;
esac
shortloc="/home/spalding/Desktop"
cp $fileloc/*.desktop $shortloc/
echo "shortcuts for each scanner is created in $shortloc"






