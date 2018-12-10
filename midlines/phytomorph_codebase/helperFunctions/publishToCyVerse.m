function [] = publishToCyVerse(pubtype)
    if nargin == 0
        functionFile = '-a decode_qr.m -a mj_getFields.m -a getPlantCell_and_Mask.m -a mergeCropBoxes.m -a getPlantBoxesForFile.m -a checkerBoardAnalysis.m -a getFieldForFileName.m -a getCheckBoardfileNames.m -a getDayfileNames.m -a parseAndLabelStack.m -a overHead_main.m -a generalizeFeatureExtractor.m -a generalizeLoader.m -a applyAllLayers.m -a trainAIlayer.m -a cRunner.m -a scanalyzerMain.m -a /usr/local/MATLAB/R2017b/toolbox/images/imdata/eSFRdefaultGrayReference.mat -a /usr/local/MATLAB/R2017b/toolbox/images/imdata/eSFRdefaultColorReference.mat -a cellForDirk.m -a generatePlant.m -a generateMetaDataTree.m -a stomataMaster.m -a fftPatch.m -a myInterp2Sampler.m -a applyNetworkToPoint.m -a generateImageDomain.m -a detectEmergence_ver2.m -a rgbAndhsvExtract.m -a arborist.m -a measureCrossOver.m -a extractPhenotypesFromOverheadCamera_ver2.m -a clusterForest.m -a clusterTree.m -a clusterNode.m -a extractPhenotypesFromOverheadCamera.m -a genQRLargeFormatSheets.m -a detectEmergence -a imsubtract.m -a rgb2hsv_fast.m -a generateImageClass.m -a getPNN_func.m -a shapeVerticalStripNozzle.m -a network.m -a func_depthStack.m -a func_resizeDepthStack.m -a constantTransitionFunction.m -a myProb.m -a my_hmm.m -a hmm_node.m -a nozzleManifold.m -a cropImages_v2.m -a func_thumbNail.m -a smartMain_v4.m -a smartMain_v2.m -a sorguhmLeafAnalysis.m -a maizeSeedling_func3.m -a maizeSeedling_func2.m -a maizeSeedling_func1.m -a smartMain.m -a singleWholeCarrotStage2.m -a petLength.m -a partialFunction.m -a condorDeploy_ver0.m -a confocalRatio.m -a isolateRoots_overStack.m -a mecka.m -a singleWholeCarrotAnalyze.m -a op0.m -a singleSeedlingImage.m';
        cdir = dir('/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/Maize/maizeSeedling/nets/');
        cdir(1:2) = [];
        for e = 1:numel(cdir)
            functionFile = ['-a ' cdir(e).name ' ' functionFile];
        end
        tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeployR2017a/';
        mkdir(tmpCompileDirectory)
        tmpCompileDirectory = '/mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/';
        mkdir(tmpCompileDirectory)

        CMD = ['mcc -d ' tmpCompileDirectory ' ' functionFile ' -a  gmdistribution.m -a cJob.m -a im2single.m -m -v DEwrapper.m'];
        eval(CMD);
        pushCMD = ['iput -f /mnt/scratch1/phytomorph_dev/Extraction/Pipelines/deDeploy/DEwrapper /iplant/home/nmiller/publicData/DEwrapper'];
        [pushR] = system(pushCMD,'-echo');
    elseif nargin == 1
        switch pubtype
            case 'gLoader'
                %%%%%%%%%%%
                % to publish generalized loader
                loaderFunc = @(X,loaderType,loaderArgs)generalizeLoader(X,loaderType,loaderArgs);
                loaderWrapper = partialFunction(loaderFunc,'generalizeLoader');
                loaderWrapper.publish();
             case 'gFeatureExtractor'
                %%%%%%%%%%%
                % to publish generalized loader
                extractFunction = @(X,extractType,extractArgs)generalizeFeatureExtractor(X,extractType,extractArgs);
                extractFunction = partialFunction(extractFunction,'generalizeFeatureExtractor');
                extractFunction.publish();
        end
    end
end





