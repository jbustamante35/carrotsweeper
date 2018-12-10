%% make figures
trSZ = 50;
%newAT = cumsum(AT,1);
%newGT = cumsum(GT==1,1);
close all
for e = 1:500
    
    % select 50
    sel = randi(size(GT,2),trSZ,1);
    sub1 = newGT(:,sel);
    sub2 = newAT(:,sel);
    
    
    [meaA(e,:) fitSigA] = measurePOP(sub1,.8);
    [meaH(e,:) fitFigH] = measurePOP(sub2,.8);  
    
    
    %{
    u1 = mean(sub1,2);
    u2 = mean(sub2,2);
    
    xlab = 1:size(u1,1);
    [J xval] = min(abs(u1 - mean(u1)));
    [para1] = fminsearch(@(X)mySigmoid_ver0(xlab',X,u1),[u1(end) J xlab(xval)]); 
    [J,yp1] = mySigmoid_ver0(xlab',para1);
    
    xlab = 1:size(u1,1);
    [J xval] = min(abs(u2 - mean(u2)));
    [para2] = fminsearch(@(X)mySigmoid_ver0(xlab',X,u2),[u2(end) J xlab(xval)]); 
    [J,yp2] = mySigmoid_ver0(xlab',para2);
    
    
    g1 = gradient(yp1);
    g2 = gradient(yp2);
    
    %{
    plot(u1,'b')
    hold on
    plot(yp1,'c')
    hold on
    plot(u2,'r')
    plot(yp2,'m')
    
   
    
    plot(g1,'k');
    plot(g2,'g');
    %}
    plot(u1,'k');
    hold on
    [ax1,l1,l2] = plotyy(xlab,yp1,xlab,g1);
    plot(u2,'k');
    hold on
    [ax2,m1,m2] = plotyy(xlab,yp2,xlab,g2);
    
    
    set(m1,'Color','g')
    waitforbuttonpress
    hold off
    drawnow
    
    %}
  
end





%%
close all
figure
scale = 60/30*100;
plot(meaH(:,2)*scale,meaA(:,2)*scale,'.')

hold on
plot(linspace(1,5,2),linspace(1,5,2),'r')
title('max percent germ per hour')


figure
scale = 2;
plot(meaH(:,3)*scale,meaA(:,3)*scale,'k.')
hold on
plot(linspace(150,250,2),linspace(150,250,2),'r')
title('delta hours to max percent germ (hours)')


figure
scale = 30/60/24;
plot((meaH(:,4)*scale),(meaA(:,4)*scale),'k.')
hold on
plot(linspace(0,2,2),linspace(0,2,2),'r')
axis([0 2 0 2])
title('duration of germination')


figure
plot(meaH(:,1)*100,meaA(:,1)*100,'k.');
hold on
plot(linspace(0,100,2),linspace(0,100,2),'r')
title('total germ')
axis([0 100 0 100])


%%
toMSCORE = {};
toMSCORE{end+1} = Z1;
toMSCORE{end+1} = gradient(toMSCORE{end});
toMSCORE{end+1} = Z3;
toMSCORE{end+1} = gradient(toMSCORE{end});
toMSCORE{end+1} = bsxfun(@times,Z1,max(Z1,[],1).^-1);
toMSCORE{end+1} = gradient(toMSCORE{end});
toMSCORE{end+1} = bsxfun(@times,Z3,max(Z3,[],1).^-1);
toMSCORE{end+1} = cumsum(toMSCORE{end});


toMG = {};

toMG{end+1} = Z1;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = Z3;%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});
toMG{end+1} = bsxfun(@times,Z1,max(Z1,[],1).^-1);%X2;%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = gradient(toMG{end});%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = bsxfun(@times,Z3,max(Z3,[],1).^-1);%bsxfun(@minus,X,mean(X,1));
toMG{end+1} = cumsum(toMG{end});%bsxfun(@minus,X,mean(X,1));


%%
trSZ = 50;
meaA = [];
meaH = [];
for e = 1:500
    
    % select 50
    sel = randi(size(toMSCORE{1},2),trSZ,1);
    for i = 1:numel(toMSCORE)
        R1{i} = toMSCORE{i}(:,sel);
    end
    
    for i = 1:numel(toMG)
        R2{i} = toMG{i}(:,sel);
    end
    
    [germ frame] = score(R1,R2,wUG,wEG,wUF,wEF,netG{midx},net,nf);
    
    
    TRTMP = FRTOT(:,sel);
    
    [meaA(e,:) fitSigA] = measurePOP(frame,.8);
    [meaH(e,:) fitFigH] = measurePOP(TRTMP,.8);  
    
end
%% 
close all
pv = [];
for trSZ = 10:5:100
    
    for e = 1:100

        % select 50
        sel = randi(size(toMSCORE{1},2),trSZ,1);
        for i = 1:numel(toMSCORE)
            R1{i} = toMSCORE{i}(:,sel);
        end

        for i = 1:numel(toMG)
            R2{i} = toMG{i}(:,sel);
        end

        [germ frame] = score(R1,R2,wUG,wEG,wUF,wEF,netG{midx},net,nf);


        TRTMP = FRTOT(:,sel);

        [meaA(e,:) fitSigA] = measurePOP(frame,.8);
        [meaH(e,:) fitFigH] = measurePOP(TRTMP,.8);  

    end
    [J,pk] = ttest2(meaA,meaH,'dim',1);
    pv = [pv ; pk];
    plot(pv)
    drawnow
end







  