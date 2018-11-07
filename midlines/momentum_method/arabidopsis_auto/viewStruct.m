function [] = viewStruct(S,fileStack)
    for e = 1:5:numel(fileStack)
        I = imread(fileStack{e});
        I = handleFLIP(I,[]);
        hold off
        imshow(I,[]);
        hold on;
        for m = 1:numel(S{e})
            plot(S{e}(m).midlines(1,:).data(1,:),S{e}(m).midlines(1,:).data(2,:),'r');
            plot(S{e}(m).contours(1,:).data(1,:),S{e}(m).contours(1,:).data(2,:),'g');
        end
        drawnow
    end
end