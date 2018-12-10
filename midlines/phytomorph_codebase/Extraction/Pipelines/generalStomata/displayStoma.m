function [] = displayStoma(patch,para,p,cP,CL,markerType,renderSize)
    if nargin <= 5
        renderSize = [2 5];
        markerType = '.';
    end
    if nargin <= 6
        renderSize = [2 5];
    end

    if ~isempty(patch)
        imshow(patch,'border','tight','DisplayRange',[]);
        hold on
    end
    % move to the point p
    para(4:5) = p;
    plotCircle(para,CL,renderSize(1));
   
    plotCenter(para,CL);
    if ~isempty(cP)
        plotPoints(cP,p,markerType,renderSize(2));
    end
end


function [] = plotCenter(para,CL)
    plot(para(4),para(5),['k*']);
end

function [] = plotCircle(para,CL,renderSize)
    % generate transformation
    [Tx] = generateDomainTransformation(para);
    % construct the circle
    TH = linspace(0,2*pi,50);
    dispX = [cos(TH)' sin(TH)' ones(50,1)];
    dispX = (Tx*dispX')';
    % plot
    plot(dispX(:,1),dispX(:,2),CL,'LineWidth',renderSize)
    plot(dispX(1,1),dispX(1,2),[CL 'o'])
end



function [] = plotPoints(cP,p,markerType,markerSize)
    CL = {'r','g','b','m'};
    if ~isempty(cP)
        for e = 1:size(cP,2)
            plot(cP(1,e)+p(1),cP(2,e)+p(2),[CL{e} markerType],'MarkerSize',markerSize);
        end
    end
end