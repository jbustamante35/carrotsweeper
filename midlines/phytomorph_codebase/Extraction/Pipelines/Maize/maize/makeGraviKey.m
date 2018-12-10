function [key] = makeGraviKey(fileName)
    % get the file parts
    [p,n,ex] = fileparts(fileName);
    sidx = strfind(fileName,filesep);
    n = fileName(sidx(end-1)+1:(end-1));
    % make lower
    %n = lower(n);
    %{
    % fix
    if strcmp(n(end-3:end),'data')
        n = n(1:end-5);
    end
    % fix
    if ~isempty(strfind(n,'spaldingdata'))
        fidx = strfind(n,'_s');
        n = n((fidx(end)+1):end);
    end
    %}


    
    n = [n '_'];
    fidx = strfind(n,'_');
    tmpPlateName = n((fidx(1)+1):(fidx(2)-1));
    tmpPlateName = makeKey(tmpPlateName);
    fprintf([fileName '-->' tmpPlateName '\n']);


    wellName = n((fidx(2)):end);
    fidx = strfind(wellName,'_');
    wn = {};
    for i = 1:(numel(fidx)-1)
        wn{i} = wellName((fidx(i)+1):(fidx(i+1)-1));
    end
    
    for e = 1:numel(wn)
        tmp_wn = upper(wn{e});
        if numel(tmp_wn) < 2
            tmp_wn = '00';
        end
        tmp_wn = tmp_wn(1:2);
        key{e} = [tmpPlateName '*_' tmp_wn(1:end-1) '_' tmp_wn(end)];
    end
end