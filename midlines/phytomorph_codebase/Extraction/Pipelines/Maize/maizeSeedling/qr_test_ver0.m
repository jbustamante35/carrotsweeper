fileName = '/home/nate/Downloads/qrcode_test2c.tif';
I = imread(fileName);
%%
J = imcrop(I);
%%
t = imcrop(J);
msg = decode_qr(t)