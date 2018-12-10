close all
for e = 1:numel(SET)
    I = imread(SET{e}{1});
    M{e} = getMidline(SET{e}{1});
    I = imread(SET{e}{1});
    imshow(I,[]);
    hold all;
    for m = 1:numel(M{e})
        path = M{e}{m};
        path = imfilter(path,fspecial('average',[100 1]),'replicate');
        plot(path(:,1),path(:,2),'LineWidth',1);
        path = arcLength(path,'arcLen');
        NUM(m) = size(path,1);
    end
    np = min(NUM);
    midlineStack = [];
    for m = 1:numel(M{e})
        path = M{e}{m};
        path = imfilter(path,fspecial('average',[100 1]),'replicate');
        path = arcLength(path,'arcLen');
        midlineStack = cat(3,midlineStack,path(1:np,:));
    end
    midlineStack = mean(midlineStack,3);
    plot(midlineStack(:,1),midlineStack(:,2),'c','LineWidth',2);
    drawnow
    hold off
end