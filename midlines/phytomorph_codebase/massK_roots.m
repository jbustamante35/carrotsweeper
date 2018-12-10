FilePath = '/mnt/tetra/nate/MATforshiheng/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
for e = 1:numel(d)
    fn = d{e};
    [pth,nm,ext] = fileparts(d{e});
    if strcmp(ext,'.mat')
        CMD = ['iget -f "' d{e} '" /mnt/tetra/nate/MATforshiheng/']
        system(CMD)
    end
end
%%
LEN = {};
close all; 
MK = [];
cnt = 1;
for f = 1:numel(FileList)
    flag(f) = true;
    a = load(FileList{f});
    numFrames = numel(a.out);
    numFrames = 400;
    numRoots
    K = [];
    
    
    
    rootCount = [];
    for fr = 1:numFrames
        if ~isfield(a.out{fr},'ME')
            rootCount(fr) =  numel(a.out{fr}.midlines);
        else
            rootCount(fr) = NaN;
        end
    end
    
    if min(rootCount) ~= max(rootCount) & ~any(isnan(rootCount))
        flag(f) = false;
    end
    
    
    
    if flag(f)
        LEN = [];
        numRoots = min(rootCount);
        for root = 1:numRoots
    
            for fr = 1:numFrames
                if ~isfield(a.out{fr},'ME')
                    tmpM = a.out{fr}.midlines(root).data;
                    [S,LEN(root,fr)] = arcLength(tmpM','arcLen');
                    %{
                    out = cwtK_imfilter(S(1:400,:),{11});
                    K = [K out.K(20:(end-20))];
                    %}
                    %{
                    %plot3(S(:,1),S(:,2),fr*ones(size(S,1),1),'r');
                    plot(S(:,1),S(:,2),'r');
                    hold on
                    plot(S(1,1),S(1,2),'g*');
                    hold off
                    axis equal
                    drawnow
                    %}
                    LEN;
                end


            end
             cnt = cnt + 1;
            hold off
        end
        
        
        
        
        LEN = [];
        numRoots = min(rootCount);
        for root = 1:numRoots
    
            for fr = 1:numFrames
                if ~isfield(a.out{fr},'ME')
                    tmpM = a.out{fr}.midlines(root).data;
                    [S,LEN(root,fr)] = arcLength(tmpM','arcLen');
                    
                    out = cwtK_imfilter(S(1:400,:),{11});
                    K = [K out.K(20:(end-20))];
                    
                    %{
                    %plot3(S(:,1),S(:,2),fr*ones(size(S,1),1),'r');
                    plot(S(:,1),S(:,2),'r');
                    hold on
                    plot(S(1,1),S(1,2),'g*');
                    hold off
                    axis equal
                    drawnow
                    %}
                    LEN;
                end


            end
            cnt = cnt + 1;
            hold off
        end
        
        
    end
    
    
    
    rootLEN(f) = min(LEN(:));
    
    
    MK = cat(3,MK,K);
    f
end
