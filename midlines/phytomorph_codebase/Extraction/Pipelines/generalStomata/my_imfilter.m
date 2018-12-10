function [] = my_imfilter()
end


%{
    X = rand(100,50)
    Xf = fft2(X);
    K = rand(10,30);
    Kf = fft2(K);
    Rf = ifft2(mtimesx(Xf,Kf));

    fft2(x)
    y = conv_fft2(Xf, Kf, 'wrap',1);
    

%}