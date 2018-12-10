function [] = displayResults(fileName,I,returnI,BoundingBoxes,MASK,SKELETON,phenoTypes,oPath,rPath,JSONdoc)
    try
        mkdir(oPath);
        [pth,nm,ext] = fileparts(fileName);

        
        
        fileList = {};

        h = image(I/255);
        axis([1 size(I,2) 1 size(I,1)]);
        axis off
        axis equal
        hold on
        for e = 1:numel(BoundingBoxes)
            rectangle('Position',BoundingBoxes{e},'EdgeColor','r');
        end
        cropBoxImage = [oPath nm '.tif'];
        fileList{end+1} = cropBoxImage;
        fprintf(['Saving image to:' cropBoxImage '\n']);
        saveas(gca,cropBoxImage);
        close all

        

        fprintf(['Algorith shows:' num2str(numel(returnI)) ':plant(s) for display \n'])
        for e = 1:numel(returnI)
            phenoTypes(e)
            if phenoTypes(e).isPlant
                dB = bwboundaries(logical(MASK{e}));
                for k = 1:3
                    tmp = returnI{e}(:,:,k);
                    tmp(tmp(:) > 1) = 1;
                    tmp(tmp(:) < 0) = 0;
                    returnI{e}(:,:,k) = tmp;
                end
                out = flattenMaskOverlay(returnI{e},logical(MASK{e}),.35,'r');
                image(out);
                hold on;axis off;axis equal;
                for b = 1:numel(dB)
                    plot(dB{b}(:,2),dB{b}(:,1),'b');
                end
                for p = 1:numel(phenoTypes(e).paths)
                    plot(phenoTypes(e).paths(p).d(:,2),phenoTypes(e).paths(p).d(:,1),'k');
                end
                plot(linspace((phenoTypes(e).CenterOfMass(2)-phenoTypes(e).plantWIDTH(1)/2),(phenoTypes(e).CenterOfMass(2)+phenoTypes(e).plantWIDTH(1)/2),2),[phenoTypes(e).CenterOfMass(1) phenoTypes(e).CenterOfMass(1)],'c');
                plot([phenoTypes(e).CenterOfMass(2) phenoTypes(e).CenterOfMass(2)],linspace((phenoTypes(e).CenterOfMass(1)-phenoTypes(e).plantHEIGHT(1)/2),(phenoTypes(e).CenterOfMass(1)+phenoTypes(e).plantHEIGHT(1)/2),2),'c');
                TH = linspace(-pi,pi,200);
                X = phenoTypes(e).StdDis(2)*cos(TH) + phenoTypes(e).CenterOfMass(2);
                Y = phenoTypes(e).StdDis(1)*sin(TH) + phenoTypes(e).CenterOfMass(1);
                plot(X,Y,'c');
                plot(phenoTypes(e).CenterOfMass(2),phenoTypes(e).CenterOfMass(1),'co');
                plot(30000*phenoTypes(e).vMassDistribution/sum(phenoTypes(e).vMassDistribution),1:size(out,1),'r');
                plot(1:size(out,2),30000*phenoTypes(e).hMassDistribution/sum(phenoTypes(e).hMassDistribution),'r');
                plot(phenoTypes(e).longestPath(2,:),phenoTypes(e).longestPath(1,:),'r');
                fileName = [oPath nm num2str(e) '.tif'];
                fileList{end+1} = fileName;
                fprintf(['Saving image to:' fileName '\n']);
                saveas(gca,fileName);
                close all
            end
        end



        %% spin-up JSON output format

        %linkTo = stripiTicket(fileName);
        %linkPath = stripiTicket(rPath);



        [jP,tN] = fileparts(fileName);
        [md] = stripOutMetaData(tN);

        f = fields(phenoTypes);
        m = fields(md);
        phDoc = [];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % port to rsml
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        plantS = {};
        for e = 1:numel(phenoTypes)
            tmpCurve = {};
            tmpFunc = {};
            if phenoTypes(e).isPlant
                for p = 1:numel(phenoTypes(e).paths)
                    tmpCurve{p} = fliplr(phenoTypes(e).paths(p).d(:,1:2));
                    tmpFunc{p}{1}.name = 'diameter';
                    tmpFunc{p}{1}.data = 2*phenoTypes(e).paths(p).d(:,4);
                    tmpFunc{p}{2}.name = 'RGBcolor';
                    tmpFunc{p}{2}.data = phenoTypes(e).paths(p).d(:,5:7);
                end
            end
            plantS{e}.curveSet = tmpCurve;
            plantS{e}.functionSets = tmpFunc;
        end
        
        docNode = com.mathworks.xml.XMLUtils.createDocument('rsml');
        docRootNode = docNode.getDocumentElement;
        docRootNode.setAttribute('xmlns:po','http://www.plantontology.org/xml-dtd/po.dtd');
        exampleDefsCell = {};
        exampleDefsCell{1,1} = 'diameter';
        exampleDefsCell{1,2} = 'double';
        exampleDefsCell{1,3} = 'pixels';
        
        exampleDefsCell{2,1} = 'RGBcolor';
        exampleDefsCell{2,2} = 'vector';
        exampleDefsCell{2,3} = 'NA';
        
        [docNode] = generateMetaDataTree(docNode,exampleDefsCell);
        [docNode] = generatePlant(docNode,plantS);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % port to rsml
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % for each plant 
        for e = 1:numel(phenoTypes)
            % stack phenotypic data
            for fl = 1:numel(f)
                phDoc(e).phenotypes.(f{fl}) = phenoTypes(e).(f{fl});
            end
            % stack meta data
            for ml = 1:numel(m)
                phDoc(e).metaData.(m{ml}) = md.(m{ml});
            end
        end
        

        JSON_string = savejson('seedlingDoc',phDoc);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save JSON string
        fileList{end+1} = [oPath nm '_jdoc.j2son'];
        fileID = fopen(fileList{end},'w');
        fprintf(fileID,strrep(JSON_string,'\/','\\/'));
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        JSON_string = savejson('seedlingDoc',JSONdoc);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save JSON string
        fileList{end+1} = [oPath nm '_jdoc.json'];
        fileID = fopen(fileList{end},'w');
        fprintf(fileID,strrep(JSON_string,'\/','\\/'));
        fclose(fileID);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save RSML string
        fileList{end+1} = [oPath nm '_RSML.rsml'];
        xmlwrite(fileList{end},docNode);
        type(fileList{end});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save RSML string
        fileList{end+1} = [oPath nm '_devData.mat'];
        save(fileList{end},'returnI','MASK');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        % save cellmat
        %fileList{end+1} = [oPath nm '_phenotypes.csv'];
        %cell2csv(fileList{end},D);

        pushToiRods(rPath,fileList);


catch ME
    getReport(ME)
end


    
end


        %{
        % convert single number phenotypes to cells and save
        D = {};
        f = fields(phenoTypes);
        f = f([1 3 4 5 6]);
        for g = 1:numel(f)
            D{1,g} = f{g};
        end
        for e = 1:numel(phenoTypes)
            for g = 1:numel(f)
                if phenoTypes(e).isPlant
                    D{e+1,g} = phenoTypes(e).(f{g});
                else
                    D{e+1,g} = NaN;
                end

            end
        end
        %}
