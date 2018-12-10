function [MS rL] = createMasterFile(mPath)
    FileList = {};
    FileExt = {'csv'};
    verbose = 1;
    FileList = gdig(mPath,FileList,FileExt,verbose);
    Q{1} = 'para';
    Q{2} = 'error';
    Q{3} = 'swell';
    Q{4} = 'fit';
    Q{5} = 'area';
    % separaate into file types
    for e = 1:numel(Q)
        [L{e} FileList] = filterFiles(FileList,Q{e});
        [sL{e}] = getShortNameList(L{e});
    end
    
    % get short names
    for e = 1:numel(sL)
        rL{e} = {};
    end
    for i = 1:numel(sL)
        for j = 1:numel(sL)
            [J,idx1,idx2] = intersect(sL{i},sL{j});
            rL{i}{end+1} = setdiff(L{i},L{i}(idx1));
            sL{i} = sL{i}(idx1);
            sL{j} = sL{j}(idx2);
            
            L{i} = L{i}(idx1);
            L{j} = L{j}(idx2);
        end
    end
    % get meta data information
    for e = 1:numel(L)
        MD{e} = extractMetaFromFileList(L{e});
        [MD{e},idx] = sortrows(MD{e});
        L{e} = L{e}(idx);
        sL{e} = sL{e}(idx);
    end
    
    MS = {};
    % find the number of files to loop over
    SZ = numel(L{1});
    for tr = 1:SZ
        % loop over each file type
        dataVEC = [];
        H = {};
        for e = 1:numel(L)
            data = csvread(L{e}{tr});
            try
                dataVEC = [dataVEC data];
            catch
                dataVEC = [dataVEC data'];
            end
        end
        dataVEC = [(1:size(dataVEC,1))' dataVEC];
        % create header with repeating meta data
        for e = 1:size(dataVEC,1)
            for h = 1:size(MD{1},2)
                H{e,h} = MD{1}{tr,h};
            end
        end
        % move data into the header var
        dataVEC = num2cell(dataVEC);
        OFFSET = size(H);
        for h1 = 1:size(dataVEC,1)
            for h2 = 1:size(dataVEC,2)
                H{h1,OFFSET(2)+h2} = dataVEC{h1,h2};
            end
        end
        OFFSET = size(MS);
        for h1 = 1:size(H,1)
            for h2 = 1:size(H,2)
                MS{OFFSET(1)+h1,h2} = H{h1,h2};
            end
        end
    end
    colH = 1;
    
    %{
    % if view
    % build up error distribution
    ERR = [];
    for tr = 1:numel(L{1})
        type0 = 2;
        data0 = csvread(L{type0}{tr});
        ERR = [ERR;data0(:)];
    end
    eyeThresh = 1*10^-3;
    for tr = 1:numel(L{1})
        type0 = 2;
        data0 = csvread(L{type0}{tr});
        fidx = abs(data0) < eyeThresh;
        type1 = 4;
        data1 = csvread(L{type1}{tr});
        type2 = 3;
        data2 = csvread(L{type2}{tr});
        plot(data1(:,fidx),'k');hold on;plot(data2(:,fidx),'r');
        hold on
        plot(data1(:,~fidx),'b');hold on;plot(data2(:,~fidx),'m');
        hold off
        axis([0 size(data1,1) 0 .2]);
        waitforbuttonpress
    end
    %}
end

function [shortNameList] = getShortNameList(FileList)
    for e = 1:numel(FileList)
        [pth,nm,ext] = fileparts(FileList{e});
        fidx = strfind(nm,'--');
        shortNameList{e} = nm(fidx(1)+2:end);
    end
end

function [subList FileList] = filterFiles(FileList,moniker)
    sidx = [];
    for e = 1:numel(FileList)
        if ~isempty(findstr(FileList{e},moniker))
            sidx = [sidx e];
        end
    end
    subList = FileList(sidx);
    FileList(sidx) = [];
end

function [meta] = extractMetaFromFile(fileName)
    [pth nm ext] = fileparts(fileName);
    fidx1 = strfind(nm,'--');
    fidx2 = strfind(nm,'_');
    dow = nm(fidx1(1)+2:fidx2(1)-1);
    mon = nm(fidx2(1)+1:fidx2(2)-1);
    dom = nm(fidx2(2)+1:fidx2(3)-1);
    if isempty(dom)
        fidx2(3) = [];
        dom = nm(fidx2(2)+1:fidx2(3)-1);
        dom(1) = '0';
    end
    hms = nm(fidx2(3)+1:fidx2(4)-1);
    yr = nm(fidx2(5)+1:fidx1(2)-1);
    scan = nm(fidx1(2)+2:fidx1(3)-1);
    row = nm(fidx1(3)+2:end);
    if numel(row) == 1;
        row = ['0' row];
    end
    meta{1} = dow;
    meta{2} = mon;
    meta{3} = dom;
    meta{4} = hms;
    meta{5} = yr;
    meta{6} = scan;
    meta{7} = row;
end

function [meta] = extractMetaFromFileList(FileList)
    meta = {};
    for e = 1:numel(FileList)
        tM = extractMetaFromFile(FileList{e});
        for c = 1:numel(tM)
            meta{e,c} = tM{c};
        end
    end
end


%{
    mPath = '/mnt/snapper/kernelSwellingData/Scott/return/';
    mPath = '/mnt/snapper/kernelSwellingData/Scott/return_seedSource/';
    mPath = '/mnt/snapper/kernelSwellingData/Het/return/';
    mPath = '/mnt/snapper/kernelSwellingData/Scott/return_second_WIDIV/';
   



    %{
    oPath = '/mnt/snapper/kernelSwellingData/Scott/';
    fileName = 'secondWIDIV.csv';
    outName = [oPath fileName];
    cell2csv(outName,MS);
    %}
%}