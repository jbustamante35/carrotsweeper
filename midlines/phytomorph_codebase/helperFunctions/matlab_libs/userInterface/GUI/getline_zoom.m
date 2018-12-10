%Tristan Ursell
%November 2012
%Active zoom for polyline selection.
%
% [X,Y]=getline_zoom(Im1);
% [X,Y]=getline_zoom(Im1,'plot');
%
% Description:  Ever wanted to carefully select points from an image but
% could not get close enough to see where you wanted to select?  Then this
% function is for you -- it allows the user to select points from an image
% at a user specified zoom level that moves along with the points selected
% in a centered frame of reference.
%
% Im1 is the input image (matrix) from which points will be selected.
%
% Hold 'shift' or 'control' and then click to end polyline selection.
%
% This script requires the 'ginput' function.
%
% This script will generate its own figure handle.
%
% SEE ALSO: getpts, getline, getrect, zoom
%
% Example:
%
% Im1=imread('rice.png');
% [X,Y] = getline_zoom(Im1,'plot');
%

function [X,Y]=getline_zoom(Im1,varargin)

%get image size
[Sy,Sx]=size(Im1);

%check image size
if min([Sy,Sx])<10
    error('Image is too small -- just maximize it, and use ''getpts''.')
end

%create a figure handle
f0=figure;
imagesc(Im1)
xlabel('X')
ylabel('Y')
axis equal tight
title('Select an initial zoom rectangle')
colormap(gray)

%choose first zoom point
rect1=getrect;

while or(rect1(3)==0,rect1(4)==0) 
    rect1=getrect;
end

X0=round(rect1(1)+rect1(3)/2);
Y0=round(rect1(2)+rect1(4)/2);
win1=round(max(rect1(3:4)));

%start active zoom
num_pts=0;

X=[];
Y=[];

while true
    %find proper zoom bounds
    xl=round(X0-win1/2);
    xl(xl<1)=1;
    
    xr=round(X0+win1/2);
    xr(xr>Sx)=Sx;
    
    yb=round(Y0-win1/2);
    yb(yb<1)=1;
    
    yt=round(Y0+win1/2);
    yt(yt>Sy)=Sy;
    
    %get temp image
    Im_temp=Im1(yb:yt,xl:xr);
    
    %show zoom plot
    hold off
    imagesc(Im_temp)
    hold on
    xlabel('X')
    ylabel('Y')
    axis equal tight
    title(['Number of points: ' num2str(num_pts) ', Choose next point.'])
    
    if num_pts>0
        %shift all points into current zoom for plotting
        X_now=X-xl+1;
        Y_now=Y-yb+1;
        
        %show last points
        plot(X_now,Y_now,'ro')
        plot(X_now,Y_now,'r-')
        drawnow
    end
    
    %select next point
    [X_temp,Y_temp,flag0]=ginput(1);
    
    %set new zoom center
    X0=X_temp+xl-1;
    Y0=Y_temp+yb-1;
    
    %save points into vector
    X=[X,X0];
    Y=[Y,Y0];
    
    %update point count
    num_pts=num_pts+1;
    
    %display warnings
    if or(or(or(X0<0,X0>Sx),Y0<0),Y0>Sy)
        warning('Point lies outside the image.')
    end
    
    %decide whether to end polyline
    if or(flag0==3,flag0==2)
        break
    end
end

%flip vectors
X=X';
Y=Y';

%final figure ('plot')
if ~isempty(varargin)
    if strcmp(varargin{1},'plot')
        hold off
        imagesc(Im1)
        hold on
        if num_pts==1
            plot(X,Y,'go')
        elseif num_pts==2
            plot(X(1),Y(1),'go')
            plot(X(end),Y(end),'ro')
            plot(X,Y,'r-')
        else
            plot(X(1),Y(1),'go')
            plot(X(end),Y(end),'ro')
            plot(X(2:end-1),Y(2:end-1),'ro')
            plot(X,Y,'r-')
        end
        xlabel('X')
        ylabel('Y')
        axis equal tight
        title(['Number of points: ' num2str(num_pts)])
    end
else
    close(f0);
end











