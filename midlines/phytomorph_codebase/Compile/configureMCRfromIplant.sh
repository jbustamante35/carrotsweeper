#!/bin/bash
################################################################################################################
echo "Create directory for MCR-v717.zip"
mkdir -p /opt/imagePhenomics/common/
echo "Start download: MCR-v717.zip"
wget -O /opt/imagePhenomics/common/v717.zip http://davos.cyverse.org/irods-rest/rest/fileContents/iplant/home/nmiller/publicData/v717.zip?ticket=HlB1NJlWX8f8mqJ
echo "Done download: MCR-v717.zip"
################################################################################################################
################################################################################################################
echo "start unzip: v717.zip"
mkdir -p /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
#unzip /opt/imagePhenomics/common/v717.zip -C /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
sudo unzip -o -q /opt/imagePhenomics/common/v717.zip -d /opt/imagePhenomics/common/MATLAB/MATLAB_Compiler_Runtime_R2012a/
echo "Done unzip: v717.zip"
