function samplePart = readPart(fidr, spatial, spectral, framesUsed, datatype)
%readPart   samplePart = readPart(fidr, spatial, spectral, framesUsed, datatype)
%Read in part of a hyperspectral datacube.
%fider = file identifier
%Read in framesUsed at a time

if datatype==12
    V = fread(fidr,spatial*spectral*framesUsed,'uint16=>uint16');
elseif datatype==2
    V = fread(fidr,spatial*spectral*framesUsed,'int16=>int16');
end
samplePart = single(reshape(V,spatial,spectral,framesUsed));
samplePart = permute(samplePart,[3,1,2]);