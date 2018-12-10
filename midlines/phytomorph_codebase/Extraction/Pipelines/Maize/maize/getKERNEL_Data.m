function [BOX_side d_BOX_front d_BOX_top MASK_side MASK_front MASK_top cIMG flags] = getKERNEL_Data(fnSIDE,fnFRONT,fnTOP)

    [pth nm1 ext] = fileparts(fnSIDE);
    [pth nm2 ext] = fileparts(fnFRONT);
    [pth nm3 ext] = fileparts(fnTOP);
    

    fidx = strfind(nm1,'_');
    TYPE1 = nm1(fidx(end)+1:end);
    
    fidx = strfind(nm2,'_');
    TYPE2 = nm2(fidx(end)+1:end);
    
    fidx = strfind(nm3,'_');
    TYPE3 = nm3(fidx(end)+1:end);
    
    tmpS = imread(fnSIDE);
    tmpF = imread(fnFRONT);
    tmpT = imread(fnTOP);
    
    
    [BOX_side MASK_side flags.STF flags.SDF] = floorFind_ver1(tmpS);
    [BOX_front MASK_front flags.FTF flags.FDF] = floorFind_ver1(tmpF);
    [BOX_top MASK_top flags.TTF flags.TDF] = floorFind_ver1(tmpT);
    
    d_BOX_front = BOX_front;
    d_BOX_front(1) = d_BOX_front(1) + size(tmpS,2);
    d_BOX_top = BOX_top;
    d_BOX_top(1) = d_BOX_top(1) + 2*size(tmpS,2);
    
    cSIDE = flattenMaskOverlay(tmpS, logical(MASK_side), .4,'g');
    cFRONT = flattenMaskOverlay(tmpF, logical(MASK_front), .4,'r');
    cTOP = flattenMaskOverlay(tmpT, logical(MASK_top), .4,'b');
    cIMG = [cSIDE cFRONT cTOP];
end
    