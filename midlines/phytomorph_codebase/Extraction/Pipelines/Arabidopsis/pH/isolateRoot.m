function [rootCurve] = isolateRoot(I)
    disp = 0;
    G = double(rgb2gray(I));
    G = imfilter(G,fspecial('gaussian',[51 51]));
    [g1 g2] = gradient(G);
    g = (g1.^2 + g2.^2).^.5;
    curve = getLevelContours(G,50);
    ridx = [curve.length] < 2000;
    curve(ridx) = [];
    
    
    E = edge(G);
    E = imclose(E,strel('disk',5));
    E = bwareaopen(E,10);
    E = bwmorph(E,'skel',inf);
    E = imclearborder(E);
    gidx = mean(G(find(E)));
    [~,sidx] = min(abs([curve.level] - gidx));
    if disp
        imshow(I,[]);
        for e = 1:numel(curve)
            plot(curve(e).data(1,:),curve(e).data(2,:))
            hold on
        end    
        plot(curve(sidx).data(1,:),curve(sidx).data(2,:),'r');
    end
    rootCurve = curve(sidx).data;
end