function [pointList] = trackMidline(midlineM,SET,WID,DST,DSS,disp)


    tSET = SET(1:DST:end);
    diskWID = round((WID)/2);
    para.pointList = fliplr(midlineM(1:DSS:end,:));
    para.domainPara{1}.type = 'disk';
    para.domainPara{1}.value{1} = [0 diskWID diskWID];
    para.domainPara{1}.value{2} = [-pi pi round(diskWID*2*pi)];
    para.THRESH =.0001;
    para.TIME = [];
    para.h = [];
    if disp
        para.h = figure;
    end
    
    [pointList] = auto_minikini(tSET,para);
end