function [boxes] = extractBoxes(wholeImage,oPath,disp)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % start cropping out carrots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        fprintf(['Starting to crop out boxes' '\n']);
        
        boxes = cropBoxes(wholeImage,oPath,disp);
        
        fprintf(['Boxes all cropped out: ' num2str(etime(clock,tm)) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % end croppiing out carrots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % start aligning boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tm = clock;
        fprintf(['Starting to align boxes' '\n']);
        
        for i = 1:numel(boxes)
             boxes{i} = alignBoxes(boxes{i},oPath,disp);
        end
        
        fprintf(['Boxes all aligned: ' num2str(etime(clock,tm)) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % end aligning carrots
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        
        
    catch ME
        getReport(ME)
    end
end
