#!/bin/bash
################################################################################################################
echo "Create directory for MCR-v717.zip"
mkdir -p /opt/imagePhenomics/common/
echo "Start download: MCR-v717.zip"
wget -O /opt/imagePhenomics/common/v717.zip http://davos.cyverse.org/irods-rest/rest/fileContents/iplant/home/nmiller/publicData/v717.zip?ticket=HlB1NJlWX8f8mqJ
echo "Done download: MCR-v717.zip"
################################################################################################################
################################################################################################################
echo "Start download: iPlant_ver0"
wget -O iPlant_ver0 https://bitbucket.org/leeatuw/repository_experiment/raw/master/iPlant_ver0
echo "Done download: iPlant_ver0"
################################################################################################################
################################################################################################################
echo "Create directory for phytoG"
mkdir -p /opt/imagePhenomics/phytoG/
echo "Start download: phytoG"
wget -O /opt/imagePhenomics/phytoG/phytoG https://bitbucket.org/leeatuw/repository_experiment/raw/master/phytoG
sudo chmod +x /opt/imagePhenomics/phytoG/phytoG
echo "Done download: phytoG"
################################################################################################################
################################################################################################################
echo "Start download: launch.sh"
wget -O launch.sh https://bitbucket.org/leeatuw/repository_experiment/raw/master/launch.sh
echo "Done download: launch.sh"
################################################################################################################
################################################################################################################
echo "Start download: run_phytoG.sh "
wget -O /opt/imagePhenomics/phytoG/run_phytoG.sh https://bitbucket.org/leeatuw/repository_experiment/raw/master/run_phytoG.sh
sudo chmod +x /opt/imagePhenomics/phytoG/run_phytoG.sh  
echo "Done download: run_phytoG.sh "
################################################################################################################
################################################################################################################
echo "start unzip: v717.zip"
mkdir -p /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
#unzip /opt/imagePhenomics/common/v717.zip -C /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
sudo unzip -o -q /opt/imagePhenomics/common/v717.zip -d /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
echo "Done unzip: v717.zip"
################################################################################################################
################################################################################################################
echo "Start export: v717.zip"
export MCR_CACHE_ROOT=$PWD/mcr_cache
echo "Done export: v717.zip"
################################################################################################################
################################################################################################################
echo "Start run: shell for compile"
mkdir -p $MCR_CACHE_ROOT
sudo chmod +x launch.sh
./launch.sh
echo "Done run: shell for compile"
