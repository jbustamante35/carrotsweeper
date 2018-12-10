#!/bin/sh
curl -H "Pragma:" --retry 30 --retry-delay 6 -o v717.zip http://proxy.chtc.wisc.edu/SQUID/ndmiller/v717.zip
curl -H "Pragma:" --retry 30 --retry-delay 6 -o icommands.x86_64.tar.bz2 http://proxy.chtc.wisc.edu/SQUID/ndmiller/icommands.x86_64.tar.bz2
unzip -q v717.zip
export MCR_CACHE_ROOT=$PWD/mcr_cache
mkdir -p $MCR_CACHE_ROOT
unset DISPLAY
tar xvfj icommands.x86_64.tar.bz2 -C $PWD
export PATH=$PATH:$PWD/icommands/
export PATH=$PATH:$PWD/
export irodsEnvFile=$PWD/.irodsEnv
export irodsAuthFileName=$PWD/.irodsA
./run_isolateRoots_overStack.sh "MATLAB_Compiler_Runtime/v717/" ""${1}"" ""${2}"" ""${3}"" ""${4}"" ""${5}"" ""${6}"" ""${7}"" ""${8}"" ""${9}""
tar cvf ${10}.tar output
