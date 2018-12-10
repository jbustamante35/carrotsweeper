function [] = singleWholeCarrotStage2(dataFolder,betaWIDTH,betaCOUNT,nn)
    
    FilePath = './output/';
    FileList = {};
    FileExt = {'csv'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    
    for e = 1:numel(FileList)
        fprintf(['Looking@ file:' FileList{e} '\n']);
        if ~isempty(strfind(FileList{e},'histogramDistanceTransform'))
            H = csvread(FileList{e});
            [pth,nm,ext] = fileparts(FileList{e});
            fidx = strfind(nm,'_');
            nm = nm(1:(fidx(end)-1));
        end
        if ~isempty(strfind(FileList{e},'topProfile'))
            P = csvread(FileList{e});
        end
    end
    
    petWIDTH = [1 H]*betaWIDTH;
    tP = P*petWIDTH^-1;
    petCOUNT= [1 tP]*betaCOUNT;
    [C] = PCA_REPROJ(P,nn.E,nn.U);
    petLENGTH = nn.func(C');
    
    petData{1,1} = 'Petiole_Width';
    petData{1,2} = 'Petiole_Count';
    petData{1,3} = 'Petiole_Length';
    petData{2,1} = petWIDTH;
    petData{2,2} = petCOUNT;
    petData{2,3} = petLENGTH;
    cell2csv(['./output/' nm '_petioleData.csv'],petData);
end