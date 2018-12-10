 D = readtext('/home/nate/Downloads/cotdataRaw_1.csv');
 %%
 data = D(4:end,2:end);
 %%
 for e = 1:size(data,1)
     for i = 1:size(data,2)
        if isempty(data{e,i})
            data{e,i} = NaN;
        end
     end
 end
 data = cell2mat(data);
 %%
 rmidx = any(isnan(data),2);
 data(rmidx,:) = [];
 
 %%
 close all
 for e = 1:size(data,2)
    osig = data(:,e);
   
    
    wholeTM = 1:size(data,1);
    frontTM = 1:30;
    frontSnip = osig(frontTM);
    frontBrob = robustfit(frontTM,frontSnip);
    
    frontFit = frontBrob(1)+frontBrob(2)*frontTM;
    
    frontNOR = frontBrob(1)+frontBrob(2)*wholeTM;
    
    frontNOR(frontTM(end)+1:end) = frontNOR(frontTM(end));
    frontNOR = imfilter(frontNOR,fspecial('average',[1 30]),'replicate');
    %{
    plot(osig);
    hold on
    plot(frontNOR,'r');
    waitforbuttonpress
    hold off
    %}
    nsig = osig - frontNOR';
    
    nsig = nsig + 5;
    
    backTM = size(data,1)-100:size(data,1);
    backSnip = nsig(backTM);
    backBrob = robustfit(backTM,backSnip);
    backFit = backBrob(1)+backBrob(2)*backTM;
    
    backNOR = backBrob(1)+backBrob(2)*wholeTM;
    backNOR(1:backTM(1)) = backNOR(backTM(1));
    backNOR = imfilter(backNOR,fspecial('average',[1 30]),'replicate');
    rm(e) = any(nsig < 0);
    %{
    plot(nsig);
    hold on
    plot(backNOR,'g');
    waitforbuttonpress
    %}
    nsig = nsig .*(backNOR').^-1;
    %nsig = nsig - 5;
    plot(nsig,'k')
    
    %plot(frontTM,frontFit,'b')
    %plot(backTM,backFit,'g')
    
    
    %msig = zeros(size(data,1),1);
    %msig(frontTM) = frontFit;
    %msig(backTM) = backFit;
    
    %nsig = osig - msig;
    
    %midSNIP(:,e) = osig(msig==0);
    
    %plot(nsig,'r')
    
    hold off
    drawnow
    %waitforbuttonpress
    ndata(:,e) = nsig;
 end
 %%
 sigF = @(L,x0,k,x)L*(1+exp(-k*(x-x0)));
 [BETA,PSI,STATS,B] = nlmefitsa(X,Y,GROUP,V,MODELFUN,BETA0)
 