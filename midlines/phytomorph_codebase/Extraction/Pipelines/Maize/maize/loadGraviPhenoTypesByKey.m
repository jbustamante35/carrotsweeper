function [vecT SZvec data] = loadGraviPhenoTypesByKey(key,ds)
    vecT = [];
    disp = 1;
    FilePath = '/mnt/spaldingdata/nate/keyDB_Feb_03_2016/';
    fileName = [FilePath key '.mat'];
    load(fileName);
    maxTM = 61;
    MP = permute(KernelSave,[1 3 2]);
    MP = reshape(MP,[size(MP,1)*size(MP,2) size(MP,3)]);
    MP = mean(MP,1);
    midlineSave = flipdim(midlineSave,2);
    
    KernelSave = KernelSave(:,:,1:maxTM);
    root_top_Save = root_top_Save(:,:,1:maxTM);
    root_bottom_Save = root_bottom_Save(:,:,1:maxTM);
    midlineSave = midlineSave(:,:,1:maxTM);
    
    for tm = 1:maxTM
        KernelSave(:,:,tm) = bsxfun(@minus,KernelSave(:,:,tm),MP);
        root_top_Save(:,:,tm) = bsxfun(@minus,root_top_Save(:,:,tm),MP);
        root_bottom_Save(:,:,tm) = bsxfun(@minus,root_bottom_Save(:,:,tm),MP);
        midlineSave(:,:,tm) = bsxfun(@minus,midlineSave(:,:,tm),MP);
    end
    if disp
        close all
        for tm = 1:maxTM
            plot(KernelSave(:,1,tm),KernelSave(:,2,tm),'r')
            hold on
            plot(root_top_Save(:,1,tm),root_top_Save(:,2,tm),'g');
            plot(root_bottom_Save(:,1,tm),root_bottom_Save(:,2,tm),'b');
            plot(midlineSave(:,1,tm),midlineSave(:,2,tm),'k');
            plot(0,0,'k*');
            hold off
            axis equal
            drawnow
            if tm == 1
                
            end
        end
    end
    %{
    data.midlineSave = midlineSave;
    data.KernelSave = KernelSave;
    data.root_top_Save = root_top_Save;
    data.root_bottom_Save = root_bottom_Save;
    data.growthRate = growthRate;
    data.tipAngle = tipAngle;
    data.K = KUR;
    %}
    
    if ~isempty(ds)
        uM = midlineSave(1:ds:end,:,:,:);
        uK = KernelSave(1:ds:end,:,:,:);
        uUP = root_top_Save(1:ds:end,:,:,:);
        uDN = root_bottom_Save(1:ds:end,:,:,:);
        uL = squeeze(growthRate);
        uA = squeeze(tipAngle);
        uKUR = KUR;

        sz1 = size(uM);
        sz2 = size(uK);
        sz3 = size(uUP);
        sz4 = size(uDN);
        sz5 = size(uL);
        sz6 = size(uA);
        sz7 = size(uKUR);

        uMr = reshape(uM,[prod(sz1(1:3)) 1]);
        uKr = reshape(uK,[prod(sz2(1:3)) 1] );
        uUPr = reshape(uUP,[prod(sz3(1:3)) 1]);
        uDNr = reshape(uDN,[prod(sz4(1:3)) 1]);
        uKURr = reshape(uKUR,[prod(sz7(1:2)) 1]);

        %vecT = [uMr' uKr' uUPr' uDNr' uL uA uKURr'];
        vecT = [uMr' uKr' uUPr' uDNr'];
        SZvec = [sz1 0 sz2 0 sz3 0 sz4];
    end
end

%{
    FilePath = '/mnt/spaldingdata/nate/keyDB_Feb_03_2016/';
    FileList = {};
    FileExt = {'mat'};
    verbose = 1;
    FileList = gdig(FilePath,FileList,FileExt,verbose);
    data = [];
    [p,n,ext] = fileparts(FileList{1});
    [data(e,:) szV] = loadGraviPhenoTypesByKey(n,50);
    data = zeros(numel(FileList),size(data,2));
    for e = 15366:numel(FileList)
        [p,n,ext] = fileparts(FileList{e});
        [data(e,:) szV] = loadGraviPhenoTypesByKey(n,50);
        wholeKey{e} = n;
        e
        numel(FileList)
    end


    uM = zeros([size(data(1).midlineSave) numel(data)]);
    uK = zeros([size(data(1).KernelSave) numel(data)]);
    uUP = zeros([size(data(1).root_top_Save) numel(data)]);
    uDN = zeros([size(data(1).root_bottom_Save) numel(data)]);
    uL = zeros([size(data(1).growthRate) numel(data)]);
    uA = zeros([size(data(1).tipAngle) numel(data)]);
    uKUR = zeros([size(data(1).K) numel(data)]);
    uK

    for e = 1:numel(data)
        uM(:,:,:,e) = data(e).midlineSave;
        uK(:,:,:,e) = data(e).KernelSave;
        uUP(:,:,:,e) = data(e).root_top_Save;
        uDN(:,:,:,e) = data(e).root_bottom_Save;
        uL(:,:,e) = data(e).growthRate;
        uA(:,:,e) = data(e).tipAngle;
        uKUR(:,:,e) = data(e).K;
        e
    end

    ds = 50;
    uM = uM(1:ds:end,:,:,:);
    uK = uK(1:ds:end,:,:,:);
    uUP = uUP(1:ds:end,:,:,:);
    uDN = uDN(1:ds:end,:,:,:);
    uL = squeeze(uL);
    uA = squeeze(uA);

    sz1 = size(uM);
    sz2 = size(uK);
    sz3 = size(uUP);
    sz4 = size(uDN);
    sz5 = size(uL);
    sz6 = size(uA);
    sz7 = size(uKUR);

    uMr = reshape(uM,[prod(sz1(1:3)) sz1(end)]);
    uKr = reshape(uK,[prod(sz2(1:3)) sz2(end)]);
    uUPr = reshape(uUP,[prod(sz3(1:3)) sz3(end)]);
    uDNr = reshape(uDN,[prod(sz4(1:3)) sz4(end)]);
    uKURr = reshape(uKUR,[prod(sz7(1:2)) sz7(end)]);
    
    vecT = [uMr;uKr;uUPr;uDNr;uL;uA;uKURr];

    rmidx = any(isnan(vecT),2);
    vecT(rmidx,:) = [];
    [S C U E L ERR LAM] = PCA_FIT_FULL(vecT,10);
   
    uC = mean(C,2);
    h1 = figure;
    h2 = figure;
    h3 = figure;
    h4 = figure;
    for e = 1:5
        for loop = 1:3
            varC = linspace(min(C(e,:)),max(C(e,:)),10);
            for f = 1:numel(varC)
                tmpC = uC;
                tmpC(e) = varC(f);
                M = PCA_BKPROJ_T(tmpC,E,U);
                str = 1;
                stp = prod(sz1(1:3));
                tmp1 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz2(1:3));
                tmp2 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz3(1:3));
                tmp3 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz4(1:3));
                tmp4 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz5(1));
                tmp5 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz6(1));
                tmp6 = M(str:stp);
                str = stp + 1;
                stp = stp + prod(sz7(1:2));
                tmp7 = M(str:stp);
                tmp1 = reshape(tmp1,sz1(1:3));
                tmp2 = reshape(tmp2,sz2(1:3));
                tmp3 = reshape(tmp3,sz3(1:3));
                tmp4 = reshape(tmp4,sz4(1:3));
                tmp7 = reshape(tmp7,sz7(1:2));



                figure(h1);
                cnt = 1;
                for tm = [1 30 61];
                    plot(tmp1(:,1,tm),-tmp1(:,2,tm),CL{cnt})
                    hold on
                    plot(tmp2(:,1,tm),-tmp2(:,2,tm),CL{cnt});
                    plot(tmp3(:,1,tm),-tmp3(:,2,tm),CL{cnt});
                    plot(tmp4(:,1,tm),-tmp4(:,2,tm),CL{cnt});
                    axis equal
                    axis([-75 300 -100 100]);
                    drawnow
                    cnt = cnt + 1;
                end

                figure(h2);
                plot(mean(gradient(tmp5))*ones(size(tmp5)));
                hold all
                axis([0 61 0 3])

                figure(h3)
                plot(tmp6*180/pi);
                hold all
                axis([0 61 -10 90]);

                figure(h4)
                mesh(interp2(tmp7(7:end,:)));
                view([0 90]);
                if loop >= 3 | f == 1
                    waitforbuttonpress
                else
                    pause(.2)
                end
                hold off
                drawnow
                
            end
            figure(h2)
            hold off
            figure(h3)
            hold off
        end
    end

    close all
    plot3(C(1,:),C(2,:),C(3,:),'.')

    toRM = [4 5 6];
    toRM_N = 50;
    RM = [];
    for e = 1:numel(toRM)
        [J,sidx] = sort(C(toRM(e),:));
        RM = [RM sidx(1:toRM_N)];
        RM = [RM sidx(end-toRM_N+1:end)];
    end
    vecT(:,RM) = [];


    

    uM = mean(uM,4);
    uK = mean(uK,4);
    uUP = mean(uUP,4);
    uDN = mean(uDN,4);
 
    


    for loop = 1:40
        for tm = 1:size(uM,3)
            plot(uK(:,1,tm),-uK(:,2,tm),'r')
            hold on
            plot(uUP(:,1,tm),-uUP(:,2,tm),'g');
            plot(uDN(:,1,tm),-uDN(:,2,tm),'b');
            plot(uM(:,1,tm),-uM(:,2,tm),'k');
            plot(0,0,'k*');
            hold off
            axis equal
            drawnow
        pause(.1)
        end
    end

    close all;
    CL = {'r' 'g' 'b'};
    cnt = 1;
    for tm = [1 30 61];
        plot(uK(:,1,tm),-uK(:,2,tm),CL{cnt})
        hold on
        plot(uUP(:,1,tm),-uUP(:,2,tm),CL{cnt});
        plot(uDN(:,1,tm),-uDN(:,2,tm),CL{cnt});
        plot(uM(:,1,tm),-uM(:,2,tm),CL{cnt});
        axis equal
        drawnow
        cnt = cnt + 1;
	end



%}