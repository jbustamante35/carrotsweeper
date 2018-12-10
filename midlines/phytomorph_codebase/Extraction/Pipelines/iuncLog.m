function [] = iuncLog()
  
    CMD = ['iquest --no-page "SELECT DATA_NAME,COLL_NAME,META_DATA_ATTR_VALUE,META_DATA_ATTR_NAME,META_DATA_ATTR_UNITS,META_DATA_MODIFY_TIME ' ...
        'WHERE META_DATA_ATTR_NAME LIKE ''ph:l:maize-cobs:1.0:stage:%'' AND COLL_NAME LIKE ''/iplant/home/gxe/maizeData/cobData%''"'];
    [o,r] = system(CMD);
    fprintf(['Done running SQL string:' num2str(toc) '\n']);
    [r] = parseRecords(r);
    
    import java.util.Map.*;
    import java.util.HashMap;
    K = HashMap();
    for e = 1:numel(r)
        key = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        dtV = K.get(key);
        key2 = r(e).META_DATA_ATTR_NAME;
        value = r(e).META_DATA_ATTR_VALUE;
        stageIDX = str2num(key2(end));
        fidx = strfind(value,':');
        
        
        
        %stageIDX = str2num(value(1:fidx(end)-1));
        stageTIME = str2num(value((fidx(end)+1):end));
        if isempty(dtV)
            dtV = NaN*ones(8,1);
        end
        dtV(stageIDX+1) = stageTIME;
        if stageIDX==0
            datetime(stageTIME,'ConvertFrom','posixtime');
        end
        
        K.put(key,dtV);
    end
    
    keySet = K.keySet;
    itr = keySet.iterator;
    cnt = 1;
    MTS = [];
    MTA = [];
    LTV = [];
    while itr.hasNext()
        key = itr.next();
        dtV = K.get(key);
        dtS = dtV(2) - dtV(1);
        dtA = dtV(3:end) - dtV(2);
        dtA = dtA / 60;
        dtS = dtS / 60;
        MTS(cnt,:) = dtS';
        MTA(cnt,:) = dtA';
        LTV(cnt,:) = dtV(1:2)';
        datetime(dtV,'ConvertFrom','posixtime');
        cnt = cnt + 1;
    end
    
    rm = any(isnan(LTV),2);
    LTV(rm,:) = [];
    figure;
    dL = sort(LTV(:,2) - min(LTV(:,1)))/60;
    
    rate = gradient(dL).^-1;
    rate = 5*gradient(imfilter(dL,fspecial('average',[20 1]),'replicate')).^-1;
    %rate = imfilter(rate,fspecial('average',[20 1]));
    plot(dL,1:numel(dL));
    figure;
    plotyy(dL,rate,dL,1:numel(dL))
    
    %figure;
    %errorbar(nanmean(MTS),nanstd(MTS));    % Plot with errorbars
    figure;
    errorbar(nanmean(MTA),(nanstd(MTA)));
   
    set(gca,'XTickLabel',{'Fist S-line','First M-line','Post-init/Pre-image-read','Post-image-read','Algorithm-Complete','Last-M-line','Last S-Line'})
    axis(gca,[0 6.5 0 6])
    
    
    
    
end