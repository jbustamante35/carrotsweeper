function [mosaic_width, mosaic_heigth, main_img]=get_mosaic_size(main_img_file)
%% Load & Present main image
main_img=imread(main_img_file);
img_hndl=imshow(main_img);
title('Main image, subject to mosaicing', 'FontSize', 14, 'Color', [0,0,0]);

prompt = {'Resize/zoom factor'};
dlg_title = 'Input for image resize function';
num_lines = 1;
def = {'1'}; %{num2str(main_img_rows),num2str(main_img_cols)};
resize_ratio_cell = (inputdlg(prompt,dlg_title,num_lines,def));
resize_ratio=str2double(resize_ratio_cell{1});

if (resize_ratio~=1)
    main_img=imresize(main_img,resize_ratio);
    img_hndl=imshow(main_img);
    title(['Resized by ',resize_ratio_cell{1},' main image']);
end
axis on;
grid on;

%% Set mosaic element size
[img_rows, img_cols, ~]=size(main_img);
mosaic_width=ceil(img_cols/14);
mosaic_heigth=ceil(img_rows/14);
mosaicSizeText=sprintf('%dX%d', mosaic_heigth, mosaic_width);
img_axis_hndl=get(img_hndl,'Parent');
choice = questdlg('Choose mosaic size',...
    sprintf('Choose by resizing the blue rectangle in image center.\n Press finish when done.'),...
    mosaicSizeText, 'Set manually', mosaicSizeText);


if strcmpi(choice,'SET MANUALLY')
   title(sprintf('Choose moisac element size. Double click when done.'),...
      'FontSize', 18, 'Color', [1, 0, 0]);
   h_rect = imrect(img_axis_hndl, [round(img_cols/2)-mosaic_heigth/2,...
       round(img_rows/2)-mosaic_width/2, mosaic_width, mosaic_heigth]);
%    setFixedAspectRatioMode(h_rect,'true'); % Dectate sqaure shape 

   addNewPositionCallback(h_rect, @textROIPosChanged); % display selection info
   wait(h_rect); %getPosition(h_rect); % [xmin ymin width height]

   pos=getPosition(h_rect); % try finding more elegant solution for this one TODO
   mosaic_width=ceil(pos(3));
   mosaic_heigth=ceil(pos(4));
   
   delete(h_rect);
%    delete(h_text);
   title('Main image, subject to mosaicing', 'FontSize', 14, 'Color', [0,0,0]);
end
axis off;
grid off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Interal servise sub-functions                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function text_ROI(pos)
% text(pos(1)+round(pos(3)),pos(2)+round(pos(4)),['Width-',num2double(3),' Height-',num2double(4),'.']);

function textROIPosChanged(pos)
pos=round(pos);
poseMsg=sprintf('Mosaic element size: [Heigh,Width]- [%dx%d]. Double click when done.',...
   pos(3), pos(4));

title(poseMsg, 'FontSize', 14, 'Color', [1,0,0]);