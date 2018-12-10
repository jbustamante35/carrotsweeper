function [toSaveContour] = renderResults(M,I,KC)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                renderResults.m 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                M:       An image to be analyzed in a matrix.
                I:       Number of cobs that are expected to be analyzed. 
                KC:      The default pixel to be considered noise relative to 1200 dpi.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        fprintf(['starting with image and results display\n']);
        % show image 
        h = image(I);
        hold on
        toSaveContour = [];
        %%%%%%%%%%%%%%%%%%%%%%%
        % render measurements to the image
        for k = 1:numel(M)            

            % plot the major and minor axis
            plot(M(k).MajorLine(:,1),M(k).MajorLine(:,2),'r');
            plot(M(k).MinorLine(:,1),M(k).MinorLine(:,2),'b');
            plot(M(k).contour(1,2),M(k).contour(1,1),'r*');
            plot(M(k).contour(:,2),M(k).contour(:,1),'c');

            % get the eth contour and project to the kernels frame        
            toSaveContour = [toSaveContour;M(k).iContour(:)'];

        end
        % title the image with the number of objects
        title(num2str(KC));     
        % format the image
        axis equal;axis off;drawnow;set(gca,'Position',[0 0 .95 .95]);
        fprintf(['ending with image and results display\n']);
    catch ME
        close all;
        getReport(ME);
        fprintf(['******error in:renderResults.m******\n']);
    end
end