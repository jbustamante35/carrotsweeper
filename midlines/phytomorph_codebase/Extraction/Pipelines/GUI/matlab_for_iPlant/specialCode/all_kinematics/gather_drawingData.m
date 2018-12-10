function [jP] = gather_drawingData(jP)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prompt = {'Sample frequency along midline (px):',...                  
                  'Disk tracking radius (px):',...
                  'Frame stack value (fr):',...
                  'Evans Percent for Growth Zone (per):',...
                  'Precent for fractional Growth (per):',...
                  'Absolute cutoff value (%/frame):'...
                  'Max strain for colormap: (%/frame):'...
                  'Max velocity for colormap (px/fr):'...
                  'Min strain for colormap: (%/frame):'...
                  'Min velocity for colormap (px/fr):'};
        dlg_title = 'Please enter kinematic values:';
        default = {'10','70','5','.70','.90','.0008','.007','6','0','0'};             % default down sample for frame
        res = inputdlg(prompt,dlg_title,1,default);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smoothing parameter for derivative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.nP = str2num(res{1});    
        para.diskRadius = str2num(res{2});
        para.width = str2num(res{3});
        para.THRESH = 10;
        para.RELcutoff = str2num(res{4});
        para.PERcutoff = str2num(res{5});
        para.ABScutoff = str2num(res{6});
        para.colorMap_strain = [str2num(res{9}) str2num(res{7})];
        para.colorMap_velocity = [str2num(res{10}) str2num(res{8})];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sample the image via user
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:numel(jP)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample the image via user
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            curve = myDraw(jP(e).inPort.fileList{1}.fileName,'message',{'please click along midline','please click on QC'});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % combine user data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % assign the midline
            iM = curve{1};
            iM = arcLength(iM,'arcLen',0);
            %%%%%%%%%%%%%%%            
            numP = round(size(iM,1)/para.nP);
            iM = arcLength(iM,'spec',numP);                
            %%%%%%%%%%%%%%%
            iM = fliplr(iM);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % store fileList
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            % attach curve to track
            jP(e).para.pointList = iM;
            % attach anchor point
            jP(e).para.anchor_point = fliplr(curve{2});
            % attach time tracking info
            jP(e).para.TIME = [];
            % distribute the common parameters
            fNames = fieldnames(para);
            for i = 1:numel(fNames)
               jP(e).para.(fNames{i}) = para.(fNames{i});
            end
            
    end
end
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prompt = {'Time downsample values (fr):',...
                  'Time window values (fr):',...
                  'Space window values (px):',...
                  'Space downsample values (px):'...
                  'Derivative values (px):'....
                  'QC window size (px):'...
                  'Disk tracking radius (px):'};
        dlg_title = 'Please enter kinematic values:';
        default = {'1','1,3,5,7','30','50,100','15,25,30,40,50','30','70'};             % default down sample for frame
        res = inputdlg(prompt,dlg_title,1,default);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time downsample value
    % number of frames to skip in the analysis
    % the "helps" or is "needed" becuase, 
    % while the tracking "needs" a high frame rate to 
    % find matching patches, BUT over "short" time intervals,
    % the tracking is noisy.
    % note: a paper NEEDS to be written, OR a paper found
    % which address biological rates, how they are intertwined
    % with the spatio-temporal resolution of the "imaging" device.
    % as well, there is a paper (maybe the same one) which needs to 
    % address a heterogenous analysis methods, for example, the method
    % should "differ" in the expansion zone VS the meristem.
    % AND should differ from trial to trial.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.frameDownsample = str2num(res{1});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % number of frames - after downsampling
    % to use in the analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.frameWindow = str2num(res{2});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smoothing parameter for derivative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.spaceWindow = str2num(res{3});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smoothing parameter for derivative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        para.space_downsample = str2num(res{4});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smoothing parameter for derivative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.derSMOOTH = str2num(res{5});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % smoothing parameter for derivative analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.anchor_point_window = str2num(res{6});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radius of disk,thresh,downsample frames,window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.diskRadius = str2num(res{7});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radius of disk,thresh,downsample frames,window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.THRESH = 10;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % radius of disk,thresh,downsample frames,window
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        para.nP = 30;
%}