#! /bin/bash
####################################
# Get Mac Address for this machine
# Get IP Address fot this machine
####################################
#
#ifconfig | grep HWaddr
# find a line with mac address and substitue spaces with + signs
# then cut at H then cut at +. The result is mac address
MAC=$(ifconfig | grep 'HWaddr' | sed -e 's/ /+/g' | cut -d 'H' -f2 | cut -d '+' -f2)

# Get IP address fot this machine
IP=$(ifconfig | grep 'inet addr' | head -1 | sed -e 's/ /+/g' | cut -d ':' -f2 |cut -d '+' -f1)

# Add '=' as a delimiter for future use of cut.
echo "Mac address for this machine is ="$MAC
echo "IP address for this machine is ="$IP

#/home/spaldings_lab/Desktop/Scripts/Scripts_main

