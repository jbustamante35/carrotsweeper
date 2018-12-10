%%
oPath = '/home/nate/Downloads/traits/';
FilePath = oPath;
FileList = {};
FileExt = {'txt'};
FileList = gdig(FilePath,FileList,FileExt,1);

%%
clear data
for e = 1:numel(FileList)
    [data{e}] = readtext(FileList{e},'\t');
    data{e} = data{e}(2:end,12:26);
end
%%
fixTraits(data(2),3);

%% stack
DATA = [];
PKs = [];

for e = 1:numel(data)
    d = data{e}(2:end,12:26);
    for i = 1:size(d,1)
        for e2 = 1:(size(d,2))
            if strcmp(d{i,e2},'NA') | isempty(d{i,e2})
                d{i,e2} = NaN;
            end
        end
    end
    
    d = cell2mat(d);
    
   
    DATA = [DATA d];
    if e == 1
        PKs(e,:) = [1 size(d,2)];
    else
        PKs(e,:) = [PKs(e-1,end)+1 PKs(e-1,end)+size(d,2)];
    end
end
%% remove NAN
rm = any(isnan(DATA),2);
bd = DATA;
%rm = any(abs(diff(d,1,2)) > 400,2);
DATA(rm,:) = [];


%% smooth

close all
hello = 1;
rm = [];
for e = 1:size(PKs,1)
    sub = DATA(:,[PKs(e,1) : PKs(e,2)]);
    sub = imfilter(sub,fspecial('average',[1 5]),'replicate');
    DATA(:,[PKs(e,1) : PKs(e,2)]) = sub;
end
DATA(any(rm,2),:) = [];
    
%% clean 
close all
hello = 1;
TH = [300 1.25*10^4 300 450 10];
rm = [];
for e = 1:size(PKs,1)
    plot(abs(diff(DATA(:,[PKs(e,1) : PKs(e,2)]),1,2))');
    rm = [rm any(abs(diff(DATA(:,[PKs(e,1) : PKs(e,2)]),1,2)) > TH(e),2)];
    waitforbuttonpress
end
DATA(any(rm,2),:) = [];
%% normalize
for e = 1:size(PKs,1)
    sub = DATA(:,[PKs(e,1) : PKs(e,2)]);
    %[sub,mu(e,:),sigma(e,:)] = zscore(sub,1,1);
    
    ssub = sort(sub,2,'descend');
    MM(e,1) = mean(min(sub,[],2));
    
    
    MM(e,2) = mean(max(sub,[],2));
    MM(e,2) = mean(ssub(:,1:5),2);
    sub = sub - MM(e,1);
    sub = sub * MM(e,2).^-1;
    %{
    sub = bsxfun(@minus,sub,MM(e,1));
    sub = bsxfun(@minus,sub,MM(e,1));
    %}
    DATA(:,[PKs(e,1) : PKs(e,2)]) = sub;
end
%%


[S C U E L ERR LAM] = PCA_FIT_FULL(DATA,3);
%%

close all
for e = 1:size(PKs,1)
    h(e) = figure;
end
for e = 1:size(bd,1)
    if sum(isnan(bd(e,:))) < 5 & sum(isnan(bd(e,:))) > 0
       rawData = bd(e,:);
       rawData = multiZbindVec(rawData,MM,PKs,1);
        
        
        [xg(e,:),fval,exitflag,output] = fminunc(@(x)myCurveDistance(x,E,U,rawData),mean(C,1));
        sim = PCA_BKPROJ(xg(e,:),E,U);



        rawData = multiZbindVec(rawData,MM,PKs,-1);
        sim = multiZbindVec(sim,MM,PKs,-1);

        for i = 1:size(PKs,1)
            rawSUB = rawData(:,[PKs(i,1) : PKs(i,2)]);
            simSUB = sim(:,[PKs(i,1) : PKs(i,2)]);
            figure(h(i));
           
             plot(rawSUB,'k')
            hold on
            plot(simSUB,'r')
            %axis([0 13 0 1500])
            title(num2str(sum(isnan(bd(i,:)))))
            waitforbuttonpress
            hold off
        end
        

       
        hold off
    end
end


%%


beta = mvregress(X,Y)
