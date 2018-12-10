#!/bin/bash
################################################################################################################
################################################################################################################
echo "Start export: MCR cacheRoot"
export MCR_CACHE_ROOT=$HOME/mcr_cache
mkdir -p $MCR_CACHE_ROOT
echo "Done export: MCR cacheRoot"
################################################################################################################
################################################################################################################
echo "Start run: shell for compile"
$HOME/phytomorph_dev/$1 /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/MATLAB_Compiler_Runtime/v717
echo "Done run: shell for compile"


