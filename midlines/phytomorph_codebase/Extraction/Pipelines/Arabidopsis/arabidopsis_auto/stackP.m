function [db curve] = stackP(fileSet,model,oPath)
    try

        init = 1;
        final = [];
        tic
        for e = 1:numel(fileSet)
            tic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the next image
            fn = fileSet{e};
            img = myReader(fn);
            img = double(img)/255;
            % remove the time stamp
            img(1,:) = [];
            % handle flip
            img = handleFLIP(img,[]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract contour and label for tip
            [curve,~,P] = extractCurves(img,model,init,final);
            % get curve
            contour = fliplr(curve(:));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract midlines
            db{e} = extractMidlines(contour,P);
            init = P{1}{2} - 60;
            final = init + 120;
            toc
            e
        end
        toc


        [img direc] = handleFLIP(img,[]);



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % root tip
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rootS = myHS_X('myHS_X');                
        % loop over all frame
        for frame = 1:numel(db)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop over all roots in the frame
            for root = 1%:numel(db{frame})
                % create root tip sequence
                if (frame == 1);rts{root} = rootTipSequence();end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % create temp root tip                        
                rootTemp = rootTip();
                % set fields for root tip
                for feature = 1:8
                    % set data
                    %data = db{frame}{root}{feature}.data;
                    %field = db{frame}{root}{feature}.field;
                    data = db{frame}{feature}.data;
                    field = db{frame}{feature}.field;
                    rootTemp.setData(field,data);    
                end
                % set the direction flag
                rootTemp.set('flipD',direc);                        
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                % file name
                %fl = stack.get(frame-1);
                fl = fileSet{frame};
                % set the image(s)
                %rootTemp.setImageName(char(fl.getFullFileName()));
                rootTemp.setImageName(fl);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % set the root tip in the sequence
                rts{root}{frame} = rootTemp;
            end
        end


        [pth nm ext] = fileparts(fileSet{1});
        fidx = strfind(pth,filesep);
        nm = pth(fidx(end)+1:end);

        save([oPath nm '.mat'],'rts');

    catch ME
        ME
        
    end
end