function [delta] = myRegMeasure(refI,refM,movI,movM,para)
    para(1) = para(1)*pi/180;
    para;
    T = [[cos(para(1)) sin(para(1)) 0];[-sin(para(1)) cos(para(1)) 0];[para(2) para(3) 1]];
    T = affine2d(T);
    T.T;
    sz = size(refI);
    


    for k = 1:size(movI,3)
        newI(:,:,k) = imwarp(movI(:,:,k),T,'OutputView',imref2d(sz(1:2)));
    end
    newM = imwarp(movM,T,'OutputView',imref2d(sz(1:2)));

    validMov = find(newM==1);
    validRef = find(refM==1);

    validIdx = intersect(validMov,validRef);


    for e = 1:size(movI,3)
        tmpRef = refI(:,:,e);
        tmpMov = newI(:,:,e);
        delta(e) = mean((tmpRef(validIdx) - tmpMov(validIdx)).^2).^.5;
    end
    delta = mean(delta)*100;
end