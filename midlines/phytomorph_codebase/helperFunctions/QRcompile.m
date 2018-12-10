function [] = QRcompile(fileList,oPath)
    doc = table;
    cnt = 1;
    for e = 1:numel(fileList)
        [p,nm,ext] = fileparts(fileList{e});
        data = readtext(fileList{e});
        data = data{1};
        fstr = strfind(data,'{');
        fstp = strfind(data,'}');
        doc.QRfileName{cnt} = nm;
        for e = 1:numel(fstr)
            key_value_pair = data((fstr(e)+1):fstp(e)-1);
            delim = strfind(key_value_pair,'_');
            key = key_value_pair(1:(delim(1)-1));
            key = strrep(key,'/','');
            key = strrep(key,'#','');
            key = strrep(key,' ','');
            key = strrep(key,'-','');
            value = key_value_pair((delim(1)+1):end);
            doc.(key){cnt} = value;
        end
        cnt = cnt + 1;
    end
    mkdir(oPath);
    fprintf(['start write table : \n']);tic
    writetable(doc,[oPath filesep date '_QRtext_compiledData.csv'])
    fprintf(['stop write table :' num2str(toc) '\n']);
end

%{
    FilePath = '/home/nate/Downloads/qrData/';
    FileList = {};
    FileExt = {'txt'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    QRcompile(FileList)
%}