function [] = quickCheck(iPlantUser,imageFile,angleThreshold,focalLengthTarget,focalLengthThreshold,fNumberTarget,fNumberThreshold,exposureTimeTarget,exposureTimeThreshold)
    nm = '';
    angleThreshold = str2num(angleThreshold);
    focalLengthTarget = str2num(focalLengthTarget);
    focalLengthThreshold = str2num(focalLengthThreshold);
    fNumberTarget = str2num(fNumberTarget);
    fNumberThreshold = str2num(fNumberThreshold);
    exposureTimeTarget = str2num(exposureTimeTarget);
    exposureTimeThreshold = str2num(exposureTimeThreshold);
    try
        basePicturePath = '/imageStore/';
        remotePicturePath = '/iplant/home/#user#/maizeData/seedlingData/';
        remotePicturePath = strrep(remotePicturePath,'#user#',iPlantUser);

        fprintf(['RemotePath:' remotePicturePath '\n']);

        if isdeployed
            CMD = 'echo $HOME';
            [status,home] = system(CMD);
            home(1) = [];
            home(end) = [];
            file1 = [filesep home filesep 'phytomorph_dev/helperFunctions/JARs/' 'core-3.2.1.jar']
            javaaddpath(file1);
            file2 = [filesep home filesep 'phytomorph_dev/helperFunctions/JARs/' 'javase-3.2.1.jar']
            javaaddpath(file2);
        end

        if isempty(imageFile)
            imageFile = [filesep home filesep '/tmpData/tmpImage.nef'];
            fprintf(['Testing on:' imageFile '\n']);
        end


        str = input('Press Enter for Picture and q for quit','s');

        while isempty(str)
            CMD = [filesep home filesep 'phytomorph_dev/Aquire/USBcameraScripts/MNsingleCaptureImage.sh'];
            fprintf(['Running Command:' CMD '\n\n']);
            [status,echoD] = system(CMD,'-echo');
            fprintf(['starting: image load, gray, and edge \n']);
            fprintf(['working on: ' imageFile '\n']);
            
            
            newName = [filesep home filesep '/tmpData/tmpImageTMP.nef'];
            CMD = ['cp "' imageFile '" "' newName '"'];
            system(CMD);
            
            % path and name
            [p nm] = fileparts(imageFile);
            info = imfinfo(imageFile);
            focalLength = info.DigitalCamera.FocalLength;
            fNumber = info.DigitalCamera.FNumber;
            exposureTime = info.DigitalCamera.ExposureTime;
            % read the image
            I = imread(newName);
            % rectifiy
            [I angle] = rectifyImage(I);
            try
                % get QR code
                nm = getQRcode(I);
            catch
                nm = [];
            end
            % print qr code to screen
            fprintf(['*******************************\n']);
            checkPass(1) = abs(angle) < angleThreshold;
            fprintf(['Angle:' num2str(angle) ':Check_Result:' num2str(checkPass(1)) ':Target' num2str(angleThreshold) '\n']);
            fprintf(['*******************************\n']);
            checkPass(2) = abs(focalLength - focalLengthTarget) <= focalLengthThreshold;
            fprintf(['focalLength:' num2str(focalLength)  ':Check_Result:' num2str(checkPass(2)) ':Target' num2str(focalLengthTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(3) = abs(fNumber - fNumberTarget) <= fNumberThreshold;
            fprintf(['fNumber:' num2str(fNumber) ':Check_Result:' num2str(checkPass(3)) ':Target' num2str(fNumberTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(4) = abs(exposureTime - exposureTimeTarget) <= exposureTimeThreshold;
            fprintf(['exposureTime:' num2str(exposureTime) ':Check_Result:' num2str(checkPass(4)) ':Target' num2str(exposureTimeTarget) '\n']);
            fprintf(['*******************************\n']);
            checkPass(5) = ~isempty(nm);
            fprintf(['QR Text:' nm '\n']);
            fprintf(['*******************************\n']);

            
            if all(checkPass)
                if ~isempty(nm)
                    newPath = [filesep home basePicturePath date filesep];
                    CMD = ['mkdir -p "' newPath '"'];
                    fprintf(['running command:' CMD '\n']);
                    system(CMD);
                    newName = [newPath nm '.nef'];
                    CMD = ['cp "' imageFile '" "' newName '"'];
                    fprintf(['running command:' CMD '\n']);
                    system(CMD);

                    tmpRemotePath = [remotePicturePath date filesep];
                    fprintf(['running command:' CMD '\n']);
                    CMD = ['imkdir -p "' tmpRemotePath '"'];
                    system(CMD);
                    remoteName =  [tmpRemotePath nm '.nef'];
                    fprintf(['running command:' CMD '\n']);
                    CMD = ['iput -f "' newName '" "' remoteName '"'];
                    system(CMD);
                end
            else
                fprintf(['FAILED checks. Please retake picture.\n']);
            end
            str = input('Press Enter for Picture and q for quit','s');
        end
    catch ME
        getReport(ME)
    end
    
    
end

%{
    testImage = '/mnt/scratch1/JUNK/plot256.nef';
    nm = quickCheck(testImage);

    mcc -m quickCheck.m -d /mnt/scratch1/phytomorph_dev/Deploy/quickCheckMaizeSeedlings/
%}