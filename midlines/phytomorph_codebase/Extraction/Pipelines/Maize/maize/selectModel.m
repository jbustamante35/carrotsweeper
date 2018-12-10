function [D] = selectModel(dataFile,model,outPort,h)
    dataS = genLoader(dataFile,0);
    delta = bsxfun(@minus,model,dataS.V);
    target = expectedKernelNumber(dataFile);
    delta = sum(delta.*delta,2);
    [Junk,sidx] = sort(delta);
    sidx = sidx(1:target);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for view points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw data for user
    close all
    I = myReader(dataS.imageStack{1}.fileName,'toGray',1);                
    imshow(I);
    hold on;
    nProps.Color = 'g';
    dataS.gD.pointList{1}.view(h,nProps);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % view graphs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for u = 1:numel(sidx)
        vProps.Color = 'r';
        dataS.G.sG{sidx(u)}.N{1}{1}.view(h,vProps);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    % spool angles
    angle = -dataS.A(sidx,:);
    [p n ext] = fileparts(dataFile);
    fn = [outPort.csvPath n '.csv'];
    [JUNK oidx] = sort(dataS.gamma(sidx,2));
    fn
    csvwrite(fn,angle);

    % spool images
    fn = [outPort.imagePath n '.tif'];
    fn
    saveas(gca,fn);

end