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
## April.8.2015.
### April.14.2015.
### date added in .log and need to relocate .log
####################################
####################################
##############################################
# read lines from given configuration text file and store information into variables
##############################################
scannerIDdevnum=$(echo $3)
#echo $scannerIDdevnum
#echo $4
truescannerID=$(echo $4)
scannerID=$(echo $5)
scannerIDroot=$(echo $6)
scannerIDport=$(echo $7)
## Modifed(newly added)
first=$(/home/spaldings_lab/Desktop/Scripts/Scripts_main/getDeviceNameForCommand.sh Perfection 1 $scannerIDdevnum)
second=$(/home/spaldings_lab/Desktop/Scripts/Scripts_main/getDeviceNameForCommand.sh Perfection 2 $scannerIDdevnum)
echo $first
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
# get the date and replace ' ' and ':' to '_' and '-'
#Network drive does not accept ':'
datename=$(echo "$(date)" | sed -e 's/ /_/g' | sed -e 's/:/-/g')
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
echo "process status: initializing log at " date > $2.log
# init scann count
##scannerCount=1
#This will read lines as long as we have some left
## Move this part to scanner_api.sh so look up for device number there
## run_scanner.sh will directly be told device number
## while read line          
## do           
##	# read a scanner line
##    	temp=$(echo $line | cut -d: -f2)
##    	pretemp=$(echo $line | cut -d: -f1)
##	# get the true scanner ID
##	truescannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2 | cut -d":" -f1)
##	# get the scanner ID
##	scannerID[$scannerCount]=$(echo $pretemp | cut -d"D" -f2)
##	# get the scanner root
##	scannerIDroot[$scannerCount]=$(echo $temp | cut -d_ -f1)
##	# get the scanner port
##	scannerIDport[$scannerCount]=$(echo $temp | cut -d_ -f2)
##	# query for the scanner device number
##	scannerIDdevnum[$scannerCount]=$(/home/spaldinglab/Desktop/production_scripts/Scripts_main/##getDeviceForPort.sh ${scannerIDroot[$scannerCount]} ${scannerIDport[$scannerCount]})
##	echo "Request scanner ID:" ${scannerID[$scannerCount]}  " at root:" ${scannerIDroot[$scannerCount]} ##":port" ${scannerIDport[$scannerCount]} ":dev" ${scannerIDdevnum[$scannerCount]}
##	((scannerCount++))
##done
##############################################
# given configureation,  create folders depending on the number of scanners 
##############################################
# `seq 1 ${#scannerIDroot[@]}`---sequence upto number of scannerIDroot 
# This is block comment I think
: <<'Block'
for index in `seq 1 ${#scannerIDroot[@]}`;
do
if [ $1 == ${truescannerID[$index]} ]; then
	j=$index
fi
done
Block
j=$1
### Move most things to scanner api and make run scanner simpler
echo "found scanner at the position:" $j
	# forder name is based on the scanner number  
	folderName[$j]=$(echo "$imagePath"/"$datename"/"$1")
	echo "Making folder:" ${folderName[$j]}
	mkdir --parents ${folderName[$j]}

# display currently running scanner and image number
for i in `seq 1 $numberofimages`;
do
	echo "Scanning at scanner number" $1 ${scannerIDroot[$j]} ":" ${scannerIDport[$j]} "@ image number" $i
	# run single scanner script
	echo "process status: starting scan " date >> $2.log
echo 
## original
#T=$(/home/spaldings_lab/Desktop/Scripts/Scripts_main/run_single_scanner.sh ${folderName[$j]}"/"$i".tif" $res ${scannerIDroot[$j]} ${scannerIDdevnum[$j]})
## modified
##### why the variables are not in right order?
T=$(/home/spaldings_lab/Desktop/Scripts/Scripts_main/run_single_scanner.sh ${folderName[$j]}"/"$i".tif" $res $first $second $scannerIDroot $scannerIDdevnum )

	echo "process status: ending scan at " date >> $2.log
	# scanning time
	Tsec=$(echo $T | cut -d "." -f1)
# compress image to new image	
# cd "${folderName[$j]}"	
# horizental differencial compress
	echo "process status: starting compression at " date >> $2.log	
	tiffcp -c lzw:2 ${folderName[$j]}"/"$i".tif" ${folderName[$j]}"/c"$i".tif"
	# remove original and leave compressed
	mv -a ${folderName[$j]}"/c"$i".tif" ${folderName[$j]}"/"$i".tif"
	echo "process status: ending compression at " date >> $2.log
	echo -ne "Scan time @scanner number" $1 ":" $T "\n"
	echo -ne "Time left to next scan for scanner number" $1 ": " $((deltaTime-Tsec)) "\n"
# I can tag name on the process but loop breaks
# exec -a "scannerID"$1 sleep $((deltaTime-Tsec))""
	echo "process status: starting sleep at " date >> $2.log
	# adjust sleep time to make the inverval regualar 
	sleep $((deltaTime-Tsec))""
	echo "process status: ending sleep at " date >> $2.log
done
