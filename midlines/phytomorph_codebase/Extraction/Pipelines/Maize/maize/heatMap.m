I = imread('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/07S-2047-01.TIF');
I(:,:,4) = [];
[D] = readtext('/mnt/spaldingdata/nate/communications/papers/paper with jeffG/B73_CV2_xycoords.csv');
header = D(1,:);
D(1,:) = [];
D(:,1:4) = [];


[c r v] = impixel(I);
X = D(:,1);
%X = (cell2mat(X)+0)*60.8;

X = cell2mat(X);
X = X - X(1);
X = X / X(end-1);
X = X * (c(2)-c(1));
X = X + c(1);

Y = D(:,2);
%Y = (cell2mat(Y)+0)*60.8;

Y = cell2mat(Y);
Y = Y - Y(1);
Y = Y / Y(end-1);
Y = Y * (r(2)-r(1));
Y = Y + r(1);



%%
K = colormap('jet');
K = colormap('winter');
%K = colormap('autumn');
%K(:,1:2) = fliplr(K(:,1:2));
K = colormap('hot');
C = cell2mat(D(:,3));
C = C - min(C);
C = C / max(C);
C = interp1(linspace(0,1,size(K,1)),K,C);

close all
imshow(I,[]);
hold on;
for e = 1:size(C,1)
    plot(X(e),Y(e),'o','MarkerSize',10,'Color',C(e,:),'MarkerEdgeColor','k','MarkerFaceColor',C(e,:));
end
colormap hot
colorbar