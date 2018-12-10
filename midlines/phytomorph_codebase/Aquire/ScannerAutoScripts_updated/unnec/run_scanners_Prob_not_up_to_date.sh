#!/bin/bash
####################################
# given configuration text file, run single scanner script 
####################################
# User input configuration 
# intput: 	resolution: X in dpi
#		deltaTime: X in seconds
#		duration: X in seconds
#		imagepath: example /home/spaldinglab/Desktop/images
#		scannerID: (RootHub in 3 digits)_(port number in a single digit)
#		scannerID example: root hub-2, port-4 ==> 002_4
####################################
####################################
##############################################
# read lines from given configuration text file and store information into variables
##############################################
# read the resolution
read line 
res=$(echo $line | cut -d: -f2)
# read the time step
read line
deltaTime=$(echo $line | cut -d: -f2)
# read the duration of scanning
read line
duration=$(echo $line | cut -d: -f2)
# calculate the number of images
numberofimages=$((duration/deltaTime))
# read the image path to save
read line
imagePath=$(echo $line | cut -d: -f2)
# get the date
datename=$(echo "$(date)" | sed -e 's/ /_/g' | sed -e 's/:/-/g')
echo $datename
echo "I got here"
# display configuration on screen
echo "scanning at resolution:" $res "dpi"
echo "scanning at deltaTime of:" $deltaTime "sec"
echo "scanning at duration of:" $duration "sec"
echo "scanning at number of images:" $numberofimages
echo "folder date name:" $datename
##############################################
# gather scanner device number and store into scannerIDdevnum - therefore read the remaining lines
##############################################
# init scann count
scannerCount=1
while read line          
do           
	# read a scanner line
    	temp=$(echo $line | cut -d: -f2)
    	pretemp=$(echo $line | cut -d: -f1)
	# get the scanner ID
	scannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2)
	# get the scanner root
	scannerIDroot[$scannerCount]=$(echo $temp | cut -d_ -f1)
	# get the scanner port
	scannerIDport[$scannerCount]=$(echo $temp | cut -d_ -f2)
	# query for the scanner device number
	scannerIDdevnum[$scannerCount]=$(/home/spaldinglab/Desktop/production_scripts/Scripts_main/getDeviceForPort.sh ${scannerIDroot[$scannerCount]} ${scannerIDport[$scannerCount]})
	echo "Request scanner ID:" ${scannerID[$scannerCount]}  " at root:" ${scannerIDroot[$scannerCount]} ":port" ${scannerIDport[$scannerCount]} ":dev" ${scannerIDdevnum[$scannerCount]}
	((scannerCount++))
done
##############################################
# given configureation,  create folders depending on the number of scanners 
##############################################
for j in `seq 1 ${#scannerIDroot[@]}`;
	do
	folderName[$j]=$(echo "$imagePath"/"$datename"/"$j")
	echo "Making folder:" ${folderName[$j]}
	mkdir --parents ${folderName[$j]}
done
# display currently running scanner and image number
for i in `seq 1 $numberofimages`;
do
	for j in `seq 1 ${#scannerIDroot[@]}`;
	do
		echo "Scanning at scanner number " $j ${scannerIDroot[$j]} ":" ${scannerIDport[$j]} "@ image number" $i
# run single scanner script
T=$(/home/spaldinglab/Desktop/production_scripts/Scripts_main/run_single_scanner.sh ${folderName[$j]}"/"$i".tif" $res ${scannerIDroot[$j]} ${scannerIDdevnum[$j]})
	Tsec=$(echo $T | cut -d "." -f1)
	# compress image to new image	
	cd "${folderName[$j]}"	
	ls -l | grep $i.tiff | cut -d " " -f5
	# horizental differencial compress
	tiffcp -c lzw:2 $i.tif $i"c".tif
	echo -ne "Scan time @scanner number" $j ":" $T "\n"
	echo -ne "Time left to next scan for scanner number" $i ": " $((deltaTime-Tsec)) "\n"
	done
	sleep $((deltaTime-Tsec))""
	#delete original 
	rm $i.tif
done

## displaying scanner ID correctly
