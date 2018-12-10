        %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over selection points and run analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:numel(dataSet)
            try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter the image
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %If = imfilter(I,fspecial('gaussian',[51 51],21),'replicate');
                %[gr1 gr2] = gradient(If);
                %gr = (gr1.^2 + gr2.^2).^.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create interpolated midline
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
               
            %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % curvilinear goeometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create curvilinear domain
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Domain = genCurvilinearDomain(iM,100,400,800);
                plotCurvilinearGrid(I,Domain,[10 10]);
                plot(iM(:,2),iM(:,1),'g');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filter the image and obtain gradient
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                If = imfilter(I,fspecial('gaussian',[51 51],21),'replicate');
                [gr1 gr2] = gradient(If);
                gr = (gr1.^2 + gr2.^2).^.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % reshape domain for sampling
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sz = size(Domain);
                rDomain = reshape(Domain,[sz(1)*sz(2) sz(3)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample raw image on curilinear geometry
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Ici = myInterp(I,rDomain);
                Ici = reshape(Ici,sz(1:2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % sample raw image on curilinear geometry
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                grci = myInterp(gr,rDomain);
                grci = bindVec(grci);
                grci = reshape(grci,sz(1:2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for viewing 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                imshow(grci);
                hold on;
                for e = 1:size(grci,1)
                    idx(e,:) = nMax(grci(e,:),100,2);
                    plot(idx(e,:),e*ones(1,2),'r.');
                    hold on;
                    drawnow
                end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % curvilinear goeometry
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
                straightKinematics_main(dataSet{e},outPath{e},iM,RAD,THRESH,SAM,derSMOOTH,WINDOW);
            catch ME
                ME
            end
        end
        %}