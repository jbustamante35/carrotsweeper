function [] = plotCurvilinearGrid(G,SAMP,I,color)
    if nargin == 3
        if ~isempty(I)
            imshow(I,[])
            hold on
        end
    end
    
    if nargin < 4
        color = 'r';
    end
    
    for e = 1:SAMP(1):size(G,1)
        plot(G(e,:,1),G(e,:,2),color);
    end
    %plot end bar
    plot(G(end,:,1),G(end,:,2),color);
    
    for e = 1:SAMP(2):size(G,2)
        plot(G(:,e,1),G(:,e,2),color);
    end
    plot(G(:,end,1),G(:,end,2),color);
end