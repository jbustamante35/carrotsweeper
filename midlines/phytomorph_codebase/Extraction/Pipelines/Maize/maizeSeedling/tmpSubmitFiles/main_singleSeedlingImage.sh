#!/bin/sh
echo "nodeArchType:"
uname -m
echo "nodeIP:"
dig +short myip.opendns.com @resolver1.opendns.com
echo "OS:"
cat /etc/lsb-release
curl -H "Pragma:" --retry 30 --retry-delay 6 -o v840.zip http://proxy.chtc.wisc.edu/SQUID/ndmiller/v840.zip
curl -H "Pragma:" --retry 30 --retry-delay 6 -o icommands.x86_64.tar.bz2 http://proxy.chtc.wisc.edu/SQUID/ndmiller/icommands.x86_64.tar.bz2
curl -H "Pragma:" --retry 30 --retry-delay 6 -o irods-icommands-4.1.9-centos-6.installer http://proxy.chtc.wisc.edu/SQUID/ndmiller/irods-icommands-4.1.9-centos-6.installer
curl -H "Pragma:" --retry 30 --retry-delay 6 -o core-3.2.1.jar http://proxy.chtc.wisc.edu/SQUID/ndmiller/core-3.2.1.jar
curl -H "Pragma:" --retry 30 --retry-delay 6 -o javase-3.2.1.jar http://proxy.chtc.wisc.edu/SQUID/ndmiller/javase-3.2.1.jar
unzip -q v840.zip
export MCR_CACHE_ROOT=$PWD/mcr_cache
mkdir -p $MCR_CACHE_ROOT
curl -H "Pragma:" --retry 30 --retry-delay 6 -o SLIBS.tar.gz http://proxy.chtc.wisc.edu/SQUID/SLIBS.tar.gz
tar -xvf SLIBS.tar.gz
LD_LIBRARY_PATH=$PWD/SS/
unset DISPLAY
chmod +x irods-icommands-4.1.9-centos-6.installer
sh irods-icommands-4.1.9-centos-6.installer $PWD/
export PATH=$PATH:$PWD/
export IRODS_PLUGINS_HOME=$PWD/icommands/plugins/
export PATH=$PATH:$PWD/icommands/
export IRODS_ENVIRONMENT_FILE=$PWD/irods_environment.json
./run_singleSeedlingImage.sh "MATLAB_Compiler_Runtime/v840/" ${1} ${2} ${3} ${4} ${5} ${6} ${7}
tar cvf ${8}.tar output
rm v840.zip
rm icommands.x86_64.tar.bz2
rm irods-icommands-4.1.9-centos-6.installer
rm core-3.2.1.jar
rm javase-3.2.1.jar
rm dcraw
rm singleSeedlingImage
rm run_singleSeedlingImage.sh
rm .irodsA
rm .irodsEnv
rm -r output
rm SLIBS.tar.gz
