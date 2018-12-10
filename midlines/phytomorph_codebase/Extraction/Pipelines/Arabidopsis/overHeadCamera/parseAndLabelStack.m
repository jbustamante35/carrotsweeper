function [typeTable,cropTableCheckerBoard,cropTableRedSquares] = parseAndLabelStack(fileList,map,tform,resEstimate,disp)
    typeTable = table;,
    cropTableCheckerBoard = table;
    cropTableRedSquares = table;
    
    try

        cntRed = 1;
        cntBlue = 1;

        for e = 1:numel(fileList)
            fprintf(['Starting file #:' num2str(e) ':' num2str(numel(fileList)) '\n']);tic
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % assign the file name
            typeTable(e,'fileName') = {fileList{e}};
            %cropTableCheckerBoard(e,'fileName') = {fileList{e}};
            %cropTableRedSquares(e,'fileName') = {fileList{e}};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % read the image
            %I = double(imread(fileList{e}))/255;
            [I] = newRecifiedImage(fileList{e},tform,resEstimate);
            % is night
            tf = isNight(I,50/255);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if is day
            if ~tf
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % assign day
                typeTable(e,'day') = {true};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get crop boxes
                %[tform,resEstimate] = checkerBoardAnalysis(I,box,disp);
                %[innerBlueSquare_Region] = getCheckerBoardData(I,map);
                [cropBOXES_Region,innerBlueSquare_Region] = getTypeCrops(I,map);
                


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % label if checkboard type
                if isempty(innerBlueSquare_Region)
                    typeTable(e,'checkerBoard') = {false};
                else
                    typeTable(e,'checkerBoard') = {true};
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the crop boxes for red squares
                for c = 1:numel(cropBOXES_Region)
                    cropTableRedSquares(cntRed,'fileName') = {fileList{e}};
                    cropTableRedSquares(cntRed,'CropBoxes') = {{[cropBOXES_Region(c).BoundingBox]}};
                    cntRed = cntRed + 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the crop boxes for blue squares
                for c = 1:numel(innerBlueSquare_Region)
                    cropTableCheckerBoard(cntBlue,'fileName') = {fileList{e}};
                    cropTableCheckerBoard(cntBlue,'CropBoxes') = {{[innerBlueSquare_Region(c).BoundingBox]}};
                    cntBlue = cntBlue + 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                if disp
                    imshow(I,[]);
                    hold on
                    if isempty(innerBlueSquare_Region)

                        for c = 1:numel(cropBOXES_Region)
                            rectangle('Position',cropBOXES_Region(c).BoundingBox,'EdgeColor','g','LineWidth',3);
                        end

                    else
                        for c = 1:numel(innerBlueSquare_Region)
                            rectangle('Position',innerBlueSquare_Region(c).BoundingBox,'EdgeColor','c','LineWidth',3);
                        end
                    end
                    hold off
                    drawnow
                end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if is night    
            else
                % assign night
                typeTable(e,'day') = {false};
            end
            
             fprintf(['Ending file #:' num2str(e) ':' num2str(numel(fileList)) ':' num2str(toc) '\n']);
            
        end
    catch ME
        getReport(ME)
    end
    

end