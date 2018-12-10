function [ERR] = graviExtractCore(matFile,m)
        ERR = 0;
         % key save location
        keySaveLoc = '/mnt/spaldingdata/nate/keyDB_Feb_03_2016/';

        load(matFile,'numberFrames','numberSeedlings','CROPBOX');
        tipAngle = [];
        stackFlag = 1;
        tmp_toKeep = [];
        fdispi = 0;
        % down sample the frames by this amount
        frDS = 1;
        % if disp is on for rendering data to screen
        dispi = 0;
        % if disp is on for the first frame
        fdispi = 0;
        fsnipK = 1;
        snipK = 40;
        try
            
            %%%%%%%%%%%%%%%%%%%%%
            % regenerate file name and extract wells
            %%%%%%%%%%%%%%%%%%%%%
            names = {};
            tmpPlateName = '';
            [pth,nm,ext] = fileparts(matFile);
            % replace mat file name with space
            imagePath = strrep(nm,'SPACE',' ');
            % replace mat file name with slash
            imagePath = strrep(imagePath,'SLASH',filesep); 
            %imagePath = strrep(imagePath,[filesep ')'],'-)');
            imagePath = [imagePath filesep];
            fprintf([imagePath '-->' num2str(m)]);
            fidx = findstr(imagePath,filesep);
            clip = imagePath(fidx(end-1)+1:end-1);
            SETUP = clip(1:2);
            clip(1:3) = [];
            clip = ['_' clip '_'];
            fidx = findstr(clip,'_');
            for e = 2:numel(fidx)-1
                tmpNM = clip(fidx(e)+1:fidx(e+1)-1);
                if numel(tmpNM) == 2
                    names{end+1} = tmpNM;
                else
                    tmpPlateName = tmpNM;
                end
            end
            tmpPlateName = clip(fidx(1)+1:fidx(2)-1);
            %%%%%%%%%%%%%%%%%%%%%
            % regenerate file name and extract wells
            %%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%
            % if disp the first frame
            %%%%%%%%%%%%%%%%%%%%%
            if fdispi 
                % load the first image
                [pth,nm,ext] = fileparts(matFile);
                imagePath = strrep(nm,'SPACE',' ');
                imagePath = strrep(imagePath,'SLASH',filesep);
                %imagePath = strrep(imagePath,'-',filesep); 
                %imagePath = strrep(imagePath,[filesep ')'],'-)');
                imagePath = [imagePath filesep];    
                iFileList = {};
                FileExt = {'TIF','tif'};
                verbose = 1;
                iFileList = gdig(imagePath,iFileList,FileExt,verbose);
                I = imread(iFileList{1});
                imshow(I,[]);
                drawnow
                waitforbuttonpress
            end
            %%%%%%%%%%%%%%%%%%%%%
            % if disp the first frame
            %%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if the number of names is empty or there is not a match for the
            % number of names then all is lost
            if numel(names) ~= numberSeedlings | isempty(names)
                return
                %allLost{end+1} = matFile;
                
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % load the information about each seedling and over each frame
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if numel(names) == numberSeedlings & ~isempty(names)
                [imagePath] = fixGraviFileName(matFile);
                
                [key] = makeGraviKey(imagePath);
                % loop over each seedling
                for s = 1:numberSeedlings
                    if ~exist([keySaveLoc strrep(key{s},'*','ASTRIST') '.mat'])
                        try
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % init the load strings for the data
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            curveLabelsString = ['seedling' num2str(s) 'curveLabels'];
                            tipAngleString = ['seedling' num2str(s) 'angle'];
                            tipIndexString = ['seedling' num2str(s) 'tipIndex'];
                            tipVecString = ['seedling' num2str(s) 'tipVec'];        
                            % call the load and pout data into fW
                            fW = load(matFile,curveLabelsString,tipAngleString,tipIndexString,tipVecString);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % find the new anchor points
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            for tm = 1:frDS:numberFrames
                                % get the labels
                                tmpLabels = fW.(curveLabelsString){tm};
                                zone2 = tmpLabels == 2;
                                zone6 = tmpLabels == 6;

                                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                                f = load(matFile,curveString);
                                f.(curveString).data = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',f.(curveString).data);

                                loc2(:,tm) = mean(f.(curveString).data(:,zone2),2);
                                loc6(:,tm) = mean(f.(curveString).data(:,zone6),2);
                            end
                            mloc2 = mean(loc2,2);
                            mloc6 = mean(loc6,2);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % load the curve for obtaining the growthrate
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tipPath = [];
                            midlineSave = [];
                            root_top_Save = [];
                            root_bottom_Save = [];
                            KernelSave = [];
                            TZ2 = [];
                            TZ6 = [];
                            KUR = [];


                            iFileList = {};
                            FileExt = {'TIF','tif'};
                            verbose = 1;
                            iFileList = gdig(imagePath,iFileList,FileExt,verbose);
                            imagePath
                            I = imread(iFileList{1});


                            for tm = 1:frDS:numberFrames

                                tmpLabels = fW.(curveLabelsString){tm};
                                zones = tmpLabels == 2 | tmpLabels == 6;

                                curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                                f = load(matFile,curveString);
                                f.(curveString).data = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',f.(curveString).data);

                                [idx2] = snapTo(f.(curveString).data',mloc2');
                                [idx6] = snapTo(f.(curveString).data',mloc6');


                                tipIDX = fW.(tipIndexString)(tm);
                                kernelC_tmp = [f.(curveString).data(:,idx2:tipIDX)];
                                kernelC_tmp_L = diff(kernelC_tmp,1,2);
                                kernelC_tmp_L = sum(kernelC_tmp_L.*kernelC_tmp_L,1).^.5;
                                root_top_Save(:,:,tm) = interp1(linspace(0,sum(kernelC_tmp_L),size(kernelC_tmp,2)),kernelC_tmp',linspace(0,(sum(kernelC_tmp_L)),5000));


                                kernelC_tmp = fliplr([f.(curveString).data(:,tipIDX:idx6)]);
                                kernelC_tmp_L = diff(kernelC_tmp,1,2);
                                kernelC_tmp_L = sum(kernelC_tmp_L.*kernelC_tmp_L,1).^.5;
                                root_bottom_Save(:,:,tm) = interp1(linspace(0,sum(kernelC_tmp_L),size(kernelC_tmp,2)),kernelC_tmp',linspace(0,(sum(kernelC_tmp_L)),5000));



                                kernelC_tmp = [f.(curveString).data(:,idx6:end) f.(curveString).data(:,1:idx2)];
                                kernelC_tmp_L = diff(kernelC_tmp,1,2);
                                kernelC_tmp_L = sum(kernelC_tmp_L.*kernelC_tmp_L,1).^.5;
                                KernelSave(:,:,tm) = interp1(linspace(0,sum(kernelC_tmp_L),size(kernelC_tmp,2)),kernelC_tmp',linspace(0,(sum(kernelC_tmp_L)),5000));




                                TZ2(:,tm) = f.(curveString).data(:,idx2);
                                TZ6(:,tm) = f.(curveString).data(:,idx6);

                                MASK = poly2mask(f.(curveString).data(1,:),f.(curveString).data(2,:),size(I,1),size(I,2));
                                skel = bwmorph(MASK,'skel',inf);
                                skelPoints = [];
                                [skelPoints(:,1),skelPoints(:,2)] = find(skel);
                                startPoint = f.(curveString).data(:,fW.(tipIndexString)(tm));
                                path = traceMaizeMidline(skelPoints,startPoint',mean([TZ2(:,tm) TZ6(:,tm)],2)');
                                path = fliplr(path);

                                tmpData = path;
                                tmpCurve = tmpData;
                                %tmpCurve = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',tmpCurve);
                                tmpData = diff(tmpData,1,2);
                                tmpData = sum(tmpData.*tmpData,1).^.5;
                                tipPath = [tipPath sum(tmpData)];
                                midlineSave(:,:,tm) = interp1(linspace(0,sum(tmpData),size(tmpCurve,2)),tmpCurve',linspace(0,ceil(sum(tmpData)),1000));
                                %{
                                [out] = cwtK_imfilter(tmpCurve',{[7]});
                                KUR(:,tm) = out.K(fsnipK:snipK);
                                %}
                            end
                            growthRate = tipPath;







                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % show the movie if dispi is on
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if dispi
                                for tm = 1:frDS:numberFrames
                                    try
                                        curveString = ['seedling' num2str(s) 'curve' num2str(tm)];
                                        f = load(matFile,curveString);
                                        f.(curveString).data = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',f.(curveString).data);
                                        midlineString = ['seedling' num2str(s) 'midline' num2str(tm)];
                                        tmpData = load(matFile,midlineString);
                                        tmpData = tmpData.(midlineString);
                                        tmpData = flipdim(tmpData,1);
                                        tmpData = bsxfun(@plus,[200 50]'+CROPBOX{s}(1:2)',tmpData);
                                        if dispi
                                            I = imread(iFileList{tm});
                                            imshow(I,[]);
                                            hold on
                                            plot(f.(curveString).data(1,:),f.(curveString).data(2,:),'w');
                                            plot(tmpData(1,:),tmpData(2,:),'g')

                                            plot(midlineSave(:,2,tm),midlineSave(:,1,tm),'r');
                                            plot(midlineSave(10,2,tm),midlineSave(10,1,tm),'m*');
                                            hold on                    
                                            UQg = unique(fW.(curveLabelsString){tm});
                                            for g = 1:numel(UQg)
                                                fidx = find(fW.(curveLabelsString){tm}==UQg(g));
                                                plot(f.(curveString).data(1,fidx),f.(curveString).data(2,fidx),[CL{UQg(g)} '*'],'MarkerSize',3);
                                            end
                                            plot(TZ2(1,tm),TZ2(2,tm),'r*');
                                            plot(TZ6(1,tm),TZ6(2,tm),'r*');

                                            plot(KernelSave(:,1,tm),KernelSave(:,2,tm),'r')
                                            plot(root_top_Save(:,1,tm),root_top_Save(:,2,tm),'k')
                                            plot(root_bottom_Save(:,1,tm),root_bottom_Save(:,2,tm),'k')
                                            hold off
                                            axis equal
                                            axis([CROPBOX{s}(1)+200 CROPBOX{s}(1)+CROPBOX{s}(3)+200 CROPBOX{s}(2)+50 CROPBOX{s}(2)+CROPBOX{s}(4)+50])
                                            drawnow
                                        end
                                        fprintf(['Done dataset ' num2str(m) ':' num2str(numel(FileList)) '\n']);
                                    catch ME
                                        stackFlag = 0;
                                    end
                                end
                            end

                            %{
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % stack the tip angle
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tipAngle = [tipAngle;fW.(tipAngleString)];
                            masterName = {masterName{:} names{s}};
                            setupName = {setupName{:} SETUP};
                            plateName = {plateName{:} tmpPlateName};
                            seedlingNumber = [seedlingNumber s];
                            matFile{end+1} = matFile;
                            iPath{end+1} = imagePath;
                            tmp_toKeep = [tmp_toKeep;1];
                            masterK = cat(3,masterK,KUR);
                            sN = [sN s];
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % stack the tip angle
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %}

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % save the data at the key
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            tipAngle = fW.(tipAngleString);
                            save([keySaveLoc strrep(key{s},'*','ASTRIST') '.mat'],'midlineSave','KUR','growthRate','tipAngle','KernelSave','root_top_Save','root_bottom_Save');


                            % catch error and display loading error
                        catch ME
                            disp(ME.getReport);
                            fprintf(['Error in ' num2str(m) '\n']);
                            ERR = 1;
                        end
                    end
                end
            end
        catch ME
            ME
            fprintf(['Error in ' num2str(m) '\n']);
            ERR = 1;
        end
end