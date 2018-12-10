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
filePath="/home/spalding/Desktop/ScannerAuto"
machType=$3
# read the resolution
read line 
res=$(echo $line | cut -d: -f2)
# read the time step
read line
deltaTime=$(echo $line | cut -d: -f2)
# read the duration of scanning
read line
duration=$(echo $line | cut -d: -f2)
#durationH=$(echo $line | cut -d: -f2)
# duration in minutes
#durationM=$((durationH/60))
# duration in seconds
#duration=$((durationM/60))
# calculate the number of images
numberofimages=$((duration/deltaTime))
# read the image path to save
read line
imagePath=$(echo $line | cut -d: -f2)
# get the date and replace ' ' and ':' to '_' and '-'
datename=$(echo "$(date)" | sed -e 's/ /_/g' | sed -e 's/:/-/g')
# read the routing information
read line
routePath=$(echo $line | cut -d: -f2)
# display configuration on screen
echo "scanning at resolution:" $res "dpi"
echo "scanning at deltaTime of:" $deltaTime "sec"
echo "scanning at duration of:" $duration "sec"
echo "scanning at number of images:" $numberofimages
echo "folder date name:" $datename
###read until line 4
##############################################
# gather scanner device number and store into scannerIDdevnum - therefore read the remaining lines
##############################################
# init log
echo "process status: initializing log" > $2.log
# init scann count
scannerCount=1
#This will read lines as long as we have some left
while read line          
do           
	# read a scanner line
    	temp=$(echo $line | cut -d: -f2)
    	pretemp=$(echo $line | cut -d: -f1)
	# get the true scanner ID
	truescannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2 | cut -d":" -f1)
	# get the scanner ID
	scannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2)
	# get the scanner root
	scannerIDroot[$scannerCount]=$(echo $temp | cut -d_ -f1)
	# get the scanner port
	scannerIDport[$scannerCount]=$(echo $temp | cut -d_ -f2)
	# query for the scanner device number
	scannerIDdevnum[$scannerCount]=$($filePath/getDeviceForPort.sh ${scannerIDroot[$scannerCount]} ${scannerIDport[$scannerCount]} $machType)
	echo "Request scanner ID:" ${scannerID[$scannerCount]}  " at root:" ${scannerIDroot[$scannerCount]} ":port" ${scannerIDport[$scannerCount]} ":dev" ${scannerIDdevnum[$scannerCount]}
	((scannerCount++))
done
##############################################
# given configureation,  create folders depending on the number of scanners 
##############################################
#`seq 1 ${#scannerIDroot[@]}`---sequence upto number of scannerIDroot 
for index in `seq 1 ${#scannerIDroot[@]}`;
do
if [ $1 == ${truescannerID[$index]} ]; then
	j=$index
fi
done
echo "found scanner at the position:" $j
	# create folder name
	folderName[$j]=$(echo "$imagePath"/"$datename"/"$1")
	# create routing folder name
	rfolderName[$j]=$(echo "$routePath"/"$datename"/"$1")
	# make folder for saving images
	echo "Making folder:" ${folderName[j]}	
	mkdir --parents ${folderName[j]}
	# make routing folder name
	mkdir --parents ${rfolderName[j]}

# display currently running scanner and image number
for i in `seq 1 $numberofimages`;
do
	echo "Scanning at scanner number" $1 ${scannerIDroot[$j]} ":" ${scannerIDport[$j]} "@ image number" $i
# run single scanner script
	echo "process status: starting scan" >> $2.log
T=$($filePath/run_single_scanner.sh ${folderName[$j]}"/"$i".tif" $res ${scannerIDroot[$j]} ${scannerIDdevnum[$j]})
	echo "process status: ending scan" >> $2.log
	Tsec=$(echo $T | cut -d "." -f1)
	# compress image to new image
	# horizental differencial compress
	echo "process status: starting compression" >> $2.log	
	tiffcp -c lzw:2 ${folderName[$j]}"/"$i".tif" ${folderName[$j]}"/c"$i".tif"
	# move the compressed image over the org image
	mv ${folderName[$j]}"/c"$i".tif" ${folderName[$j]}"/"$i".tif"
	# move the image to the routing folder
	mv ${folderName[$j]}"/"$i".tif" ${rfolderName[$j]}"/"$i".tif" 
	echo "process status: ending compression" >> $2.log
	echo -ne "Scan time @scanner number" $1 ":" $T "\n"
	echo -ne "Time left to next scan for scanner number" $1 ": " $((deltaTime-Tsec)) "\n"
	# I can tag name on the process but loop breaks
	#exec -a "scannerID"$1 sleep $((deltaTime-Tsec))""
	echo "process status: starting sleep" >> $2.log
	sleep $((deltaTime-Tsec))""
	echo "process status: ending sleep" >> $2.log
	# copy the file to the routepath
done

