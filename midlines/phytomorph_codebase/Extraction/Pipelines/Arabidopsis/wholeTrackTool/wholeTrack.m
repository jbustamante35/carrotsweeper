function [pointListL pointListE SE SL] = wholeTrack(fileList,pointList,disp,miniSKIP)
    %[pointList1 pointList2] = ndgrid(pointList(1,1):miniSKIP:(pointList(1,2)-1),pointList(2,1):miniSKIP:(pointList(2,2)-1));
    %para.pointList = [pointList2(:) pointList1(:)];
    para.pointList = pointList;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate vars
    RAD = 70;
    para.patchSize = RAD;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate domain
    domainPara{1}.type = 'disk';
    domainPara{1}.value{1} = [0 RAD RAD+1];
    domainPara{1}.value{2} = [-pi pi 200];
    para.domainPara = domainPara;
    if disp
        para.h = figure;
    else
        para.h = [];
    end
    para.TIME = [];
    para.THRESH = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform tracking
    [pointListL pointListE SE SL] = minikiniWHOLE(fileList,para);
    close all 
end

