%% load NAM from database
p='/home/nate/postgresql-9.2-1003.jdbc4.jar';
javaaddpath(p);
conn = database('maize','maizeuser','PVsgC4jmtyS6','org.postgresql.Driver', 'jdbc:postgresql://128.104.98.78/maize');
q = ['SELECT * FROM population_lines ' ...
                'JOIN kernel_plates ' ...
                'ON population_lines.id=kernel_plates.population_line_id ' ... 
                'JOIN kernels ' ...
                'ON kernel_plates.id=kernels.plate_id ' ...                
                'JOIN kernel_3d ' ...
                'ON kernels.id=kernel_3d.kernel_id ' ...             
                'JOIN reporting.kernel_dims_report_tbl '...
                'on kernels.id=reporting.kernel_dims_report_tbl.kernel_id '...
                'JOIN files ' ...
                'on kernels.img_gray_side_fid = files.id '];
cursor = exec(conn, q);
cursor = fetch(cursor);
fieldString = columnnames(cursor,1);
results = cursor.Data;
%{
img_bw_front_fid    | integer          | 
 img_bw_side_fid     | integer          | 
 img_bw_top_fid      | integer          | 
 img_gray_front_fid  | integer          | 
 img_gray_side_fid   | integer          | 
 img_gray_top_fid    | integer          | 
 img_color_front_fid
%}
%%
for e = 1:size(results,1)
    data{e} = char(results{e,3});
end
unique(data)
%%
ext = '_S2.tiff';
%ext = '_F.tiff';
BAD = {};
for e = 1:1:size(results,1)
    I = [];
    try
        % normal file name
        fileName = results{e,end};
        I = imread(fileName);
    catch
        % try adding _S2.tiff
        try
            fileName = char(results{e,end});
            I = imread([fileName ext]);
        catch
            try
                % replaing /_ with /            
                fileName = char(results{e,end});
                fileName = strrep(fileName,'/_','/');
                I = imread([fileName ext]);
            catch
                
                BAD{end+1} = results{e,end};
            end
        end
        
    end
    e
    size(results,1)
    [uB(e)] = backGroundKernel(I);
    imshow(I,[])
    drawnow
end
%%
close all
fidx = find(uB < 100);
ridx = find(uB==-1);
fidx = setdiff(fidx,ridx);
for e = 1:numel(fidx)
    fileName = results{fidx(e),end};
    
    try
        % normal file name
        fileName = results{fidx(e),end};
        I = imread(fileName);
    catch
        % try adding _S2.tiff
        try
            fileName = char(results{fidx(e),end});
            I = imread([fileName ext]);
        catch
            try
                % replaing /_ with /            
                fileName = char(results{fidx(e),end});
                fileName = strrep(fileName,'/_','/');
                I = imread([fileName ext]);
            catch                
                
            end
        end
        
    end
    I = double(I);
    E = edge(I);
    f = fspecial('gaussian',[21 21],4);
    [d1 d2] = gradient(imfilter(I,f,'replicate'));    
    S = std((I),1,1);
    dS = gradient(S);
    S = sum(E,1);
    dS = mean(d1,1);
    dS = abs(dS);
    PADA = 20;
    ddS = abs(imfilter(gradient(dS),f));
    cutOff = mean([ddS(1:50) ddS(end-50:end)]);
    BLOCK = ddS > 10*cutOff;
    R = regionprops(BLOCK,'PixelIdxList');
    BLOCK = zeros(size(BLOCK));
    BLOCK(R(1).PixelIdxList(1)+PADA) = 1;
    BLOCK(R(end).PixelIdxList(end)-PADA) = 1;    
    bkidx = find(BLOCK);
    BKI = [I(:,1:bkidx(1)-50) I(:,bkidx(2)+50:end)];
    KI = [I(:,bkidx(1):bkidx(end))];
    %I = [I BKI KI];
   
    sBK = std(BKI,1,2);
    kBK = std(KI,1,2);
    sBK = mean(BKI,2);
    sBK = mean(KI,2);
    bI = I - repmat(sBK,[1 size(I,2)]);
    SIG = std(bI,1,2);
    uSIG = mean(bI,2);
    ug1 = mean(abs(d1),2);
    ug1 = bindVec(ug1);
    SIG = bindVec(SIG);
    uSIG = bindVec(uSIG);
    SIG = SIG.*uSIG.*ug1;
    SIG = bindVec(SIG);
    %{
    gSIG = gradient(imfilter(SIG,f,'replicate'));
    gSIG = abs(gSIG);
    gSIG = bindVec(gSIG);
    SIG = gSIG;
    %}
    THRESH = mean([SIG(1:100);SIG(end-100:end)]);
    THRESH = graythresh(SIG);
    %SUG = imfilter(SIG,f,'replicate') > THRESH;
    PER = -.7;
    SUG = SIG > (THRESH + PER*THRESH);
    
    
    
    R = regionprops(SUG,'PixelIdxList','Area');
    A = [R.Area];
    [ma sidx] = max(A);
    
    %SUG = zeros(size(SUG));
    %SUG(R(sidx).PixelIdxList) = 1;
    
    suidx = find(SUG);
    HI = zeros(size(SUG));
    HI(suidx(1)) = 1;
    HI(suidx(end)) = 1;
    
    imshow(I,[])
    hold on
    %plot(abs(sBK)*10+500,1:numel(sBK),'r')
    %plot(abs(kBK)*10+500,1:numel(sBK),'b')
    fBK = abs(kBK - sBK);
    %plot(100*(fBK > 5)+500,1:numel(sBK),'g')
    
    gidx = find(fBK > 5);
    %plot(1:size(I,2),gidx(end)*ones(1,size(I,2)),'m');
    plot(abs(SIG)*100+50,1:numel(SIG),'r')
    plot(abs(SUG)*100+50,1:numel(SUG),'m')
    plot(abs(HI)*size(I,2),1:numel(SUG),'m')
    
    %{
    plot(S*10)
    plot(dS*300+100,'r')
    plot(ddS*5000+100,'g')
    %}
    plot(BLOCK*5000+100,'c')
    hold off
    drawnow
    pause(.1)    
    %{
    
    BK(:,e) = mean(I,2);
    e
    numel(fidx)
    drawnow
    %}
end
%%
close all
plot(BK)
