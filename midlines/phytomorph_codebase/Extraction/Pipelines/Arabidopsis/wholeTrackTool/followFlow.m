function [] = followFlow(I,f,pt)
    for e = 1:(size(f,4)-1)
        pt(e+1,:) = [ba_interp2(squeeze(f(:,:,1,e+1)),pt(1,2),pt(1,1)) ba_interp2(squeeze(f(:,:,2,e+1)),pt(1,2),pt(1,1))];
    end
    figure;
    imshow(I,[]);
    hold on;
    plot(pt(:,2),pt(:,1),'r')
end