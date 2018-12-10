function [rI,angle] = rectifyImage(I)
    oI = I;
    BUMP = 50;
    [msg,qrCropBox] = getQRcode(I);
    TOPloc = qrCropBox(2) + qrCropBox(4);
    TOPloc = TOPloc + BUMP;

    I(1:TOPloc,:,:)= [];

    floorCHOP = 30;
    I((end-floorCHOP):end,:,:) = [];

    % make gray scale
    G = rgb2gray(I);
    % filter the image
    G = imfilter(G,fspecial('gaussian',[31 31],11),'replicate');
    % find edge
    E = edge(G);
    % find the 90-plumb
    [H,T,R] = hough(E','Theta',linspace(-9,9,300));
    P = houghpeaks(H,3);
    lines = houghlines(E',T,R,P,'FillGap',size(E,2)/2,'MinLength',600);
    %{
        % for display
        imshow(I,[]);
        hold on
        CL = {'r' 'g' 'b'}
        for e = 1:numel(lines)
            plot([lines(e).point1(2),lines(e).point2(2)],[lines(e).point1(1),lines(e).point2(1)],CL{e})
        end
    %}
    A = [];
    for e = 1:numel(lines)
        xy = [lines(e).point1; lines(e).point2];
        xy = diff(xy,1,1);
        angle = atan2(xy(1),xy(2))*180/pi;
        A(e) = angle;
    end
    angle = mean(A);
    rI = imrotate(oI,angle,'bicubic','crop');
    % trim off X pixel from rotation
    rI(:,1:30,:) = [];
    rI(:,end-70:end,:) = [];
    rI((end-50):end,:,:) = [];
    %{
    mask = sum(rI,3)==0;
    for e = 1:4
        fidx = find(sum(mask(:,1:end/2),1) > 100);
        rI(:,1:fidx(end),:) = [];
        rI = imrotate(rI,90);
        mask = imrotate(mask,90);
    end
    %}
end