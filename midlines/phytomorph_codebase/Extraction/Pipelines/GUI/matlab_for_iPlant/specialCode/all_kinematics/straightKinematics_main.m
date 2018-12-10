function [out] = straightKinematics_main(jP)
    out = [];
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % fileList := list of files to process
        % outPath : = directory to save data to
        % pointList := Ponits to track
        % diskRadius := radius of the disk    
        % THRESH := the distance from the edge to stop tracking
        % disp := display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 0: step-up vars
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            para = jP.para;
            fileList = jP.inPort.fileList;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRACKING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % step 1: track points for each "disk" or domain
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for diskRad = 1:numel(para.diskRadius)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate vars
                RAD = para.diskRadius(diskRad);
                para.patchSize = RAD;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % generate domain
                domainPara{1}.type = 'disk';
                domainPara{1}.value{1} = [0 RAD RAD+1];
                domainPara{1}.value{2} = [-pi pi 200];
                para.domainPara = domainPara;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % perform tracking
                para.pointList = minikini(fileList,para); 
                fprintf(['done with tracking \n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % populate the root model
                rts = rootTipSequence();
                for frame = 1:(size(para.pointList,3)-1)
                    % create temp root tip
                    trt = rootTip();
                    trt.setData('midline_1',para.pointList(:,:,frame));
                    trt.midline_1.arcLengthPara();                   
                    rts{frame} = trt;
                end
                % set the flow field for the stack to the root tip sequence
                rts.setData('position_field_midline',permute(para.pointList,[1 3 2]));
                rts.saveFlow(jP.outPort);   
                profile = rts.create_flowProfile(jP.para);               
                fileList = rootTipSequence.save_flowProfile(profile,jP.outPort);
                fprintf(['starting save \n'])
                rts.saveFlow(jP.outPort);
                close all;
                %{
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % save matlab file
                    matPath = [para.outPath filesep 'matlab_data' filesep];
                    mkdir(matPath);
                    fl = [matPath '{[patchSize]_[' num2str(RAD) ']}' 'data.mat'];
                    save(fl);
                %}
            end
        catch ME
            ME
            ME.stack.line
            ME.stack.file
        end
end


