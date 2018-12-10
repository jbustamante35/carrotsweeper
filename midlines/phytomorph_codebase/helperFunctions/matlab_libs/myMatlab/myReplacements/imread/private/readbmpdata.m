function X = readbmpdata(info)
%READBMPDATA Read bitmap data
%   X = readbmpdata(INFO) reads image data from a BMP file.  INFO is a
%   structure returned by IMBMPINFO. X is a uint8 array that is 2-D for
%   1-bit, 4-bit, and 8-bit image data.  X is M-by-N-by-3 for 24-bit and
%   32-bit image data.  

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.5 $  $Date: 2011/05/17 02:27:57 $

offset = info.ImageDataOffset;
width = info.Width;
height = info.Height;
filename = info.Filename;

switch info.CompressionType

 case 'none'

  switch info.BitDepth
   case 1
    X = logical(bmpReadData1(filename, offset, width, height));
    
   case 4
    X = bmpReadData4(filename, offset, width, height);
    
   case 8
    X = bmpReadData8(filename, offset, width, height);
    
   case 16
    X = bmpReadData16(filename, offset, width, height);
    
   case 24
    X = bmpReadData24(filename, offset, width, height);
    
   case 32
    X = bmpReadData32(filename, offset, width, height);
    
  end
  
 case '8-bit RLE'
  X = bmpReadData8RLE(filename, offset, width, height);
  
 case '4-bit RLE'
  X = bmpReadData4RLE(filename, offset, width, height);
  
 case 'bitfields'
  error(message('MATLAB:imagesci:readbmpdata:bitfieldsCompression'));
  
 case 'Huffman 1D'
  error(message('MATLAB:imagesci:readbmpdata:huffman1dCompression'));
  
 case '24-bit RLE'
  error(message('MATLAB:imagesci:readbmpdata:rle24BitCompression'));
  
end
        
%%%
%%% bmpReadData8 --- read 8-bit bitmap data
%%%
function X = bmpReadData8(fname, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
paddedWidth = 4*ceil(width/4);

fid = fopen(fname,'r','ieee-le');
status = fseek(fid,offset,'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badEightBitDataOffset'));
end
X = fread(fid,paddedWidth*abs(height),'*uint8');
fclose(fid);
  
count = length(X);
if (count ~= paddedWidth*abs(height))
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    % Fill in the missing values with zeros.
    X(paddedWidth*abs(height)) = 0;
end

if height>=0
  X = rot90(reshape(X, paddedWidth, height));
else
  X = reshape(X, paddedWidth, abs(height))';
end

if (paddedWidth ~= width)
    X = X(:,1:width);
end



%%%
%%% bmpReadData8RLE --- read 8-bit RLE-compressed bitmap data
%%%
function X = bmpReadData8RLE(fname, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
paddedWidth = 4*ceil(width/4);

fid = fopen(fname,'r','ieee-le');
status = fseek(fid,offset,'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badEightBitRleDataOffset'));
end
inBuffer = fread(fid,'*uint8');
fclose(fid);

X = bmpdrle(inBuffer, paddedWidth, abs(height), 'rle8');

if height>=0
  X = rot90(X);
else
  X = X';
end
if (paddedWidth ~= width)
    X = X(:,1:width);
end



%%%
%%% bmpReadData4 --- read 4-bit bitmap data
%%%
function X = bmpReadData4(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
paddedWidth = 8*ceil(width/8);
numBytes = paddedWidth * abs(height) / 2; % evenly divides because of padding

fid = fopen(filename,'r','ieee-le');
status = fseek(fid,offset,'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badFourBitDataOffset'));
end
XX = fread(fid,numBytes,'*uint8');
fclose(fid);

count = length(XX);
if (count ~= numBytes)
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    % Fill in the missing values with zeros.
    X(numBytes) = 0;
end
XX = reshape(XX, paddedWidth / 2, abs(height));

X = repmat(uint8(0), paddedWidth, abs(height));
X(1:2:end,:) = bitslice(XX,5,8);
X(2:2:end,:) = bitslice(XX,1,4);

if height>=0
  X = rot90(X);
else
  X = X';
end
if (paddedWidth ~= width)
    X = X(:,1:width);
end


%%%
%%% bmpReadData4RLE --- read 4-bit RLE-compressed bitmap data
%%%
function X = bmpReadData4RLE(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
paddedWidth = 8*ceil(width/8);

fid = fopen(filename,'r','ieee-le');
status = fseek(fid,offset,'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badFourBitRleDataOffset'));
end
inBuffer = fread(fid,'*uint8');
fclose(fid);

if height>=0
  X = rot90(bmpdrle(inBuffer, paddedWidth, abs(height), 'rle4'));
else
  X = bmpdrle(inBuffer, paddedWidth, abs(height), 'rle4')';
end
if (paddedWidth ~= width)
    X = X(:,1:width);
end


%%%
%%% bmpReadData1 --- read 1-bit bitmap data
%%%
function X = bmpReadData1(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
paddedWidth = 32*ceil(width/32);
numPixels = paddedWidth * abs(height);  % evenly divides because of padding

% 1-bit BMP data has big-endian byte ordering
fid = fopen(filename,'r','ieee-be');
status = fseek(fid,offset,'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badOneBitDataOffset'));
end
[X, count] = fread(fid,paddedWidth*abs(height),'*ubit1');
fclose(fid);

if (count ~= numPixels)
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    % Fill in the missing values with zeros.
    X(numPixels) = 0;
end
X = reshape(X, paddedWidth, abs(height));

if height>=0
  X = rot90(X);
else
  X = X';
end

if (paddedWidth ~= width)
    X = X(:,1:width);
end

%%%
%%% bmpReadData16 --- read 16-bit bitmap data
%%%
function RGB = bmpReadData16(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
scanlineLength = 2 * ceil(width/2);
numSamples = scanlineLength * abs(height);

fid = fopen(filename,'r','ieee-le');

status = fseek(fid,offset,0);
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badSixteenBitDataOffset'));
end

X = fread(fid, numSamples, 'uint16=>uint16');
fclose(fid);

% Fill in any missing values with zeros.
if (numel(X) < numSamples)
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    X(numSamples) = 0;
end

if (height >= 0)
  X = rot90(reshape(X, scanlineLength, abs(height)));
else
  X = reshape(X, scanlineLength, abs(height))';
end

if (scanlineLength ~= width)
    X = X(:, 1:width);
end

RGB(1:abs(height), 1:width, 1) = uint8(bitslice(X,11,15));
RGB(:,:,2) = uint8(bitslice(X,6,10));
RGB(:,:,3) = uint8(bitslice(X,1,5));

%Scale data for display
RGB = bitor(bitshift(RGB,3),bitshift(RGB,-2));

%%%
%%% bmpReadData24 --- read 24-bit bitmap data
%%%
function RGB = bmpReadData24(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
scanlineLength = 4 * ceil((3 * width) / 4);
numSamples = scanlineLength * abs(height);

fid = fopen(filename, 'r', 'ieee-le');

status = fseek(fid, offset, 'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:badTwentyFourBitDataOffset'));
end

X = fread(fid, numSamples, 'uint8=>uint8');
fclose(fid);

% Fill in any missing values with zeros.
if (numel(X) < numSamples)
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    X(numSamples) = 0;
end

if (height >= 0)
  X = rot90(reshape(X, scanlineLength, abs(height)));
else
  X = reshape(X, scanlineLength, abs(height))';
end

if (width ~= scanlineLength/3)
    X = X(:, 1:(3 * width));
end

RGB(1:abs(height), 1:width, 3) = X(:,1:3:end);
RGB(:, :, 2) = X(:,2:3:end);
RGB(:, :, 1) = X(:,3:3:end);


%%%
%%% bmpReadData32 --- read 32-bit bitmap data
%%%
function RGB = bmpReadData32(filename, offset, width, height)

% NOTE: BMP files are stored so that scanlines use a multiple of 4 bytes.
scanlineLength = 4 * width;
numSamples = scanlineLength * abs(height);

fid = fopen(filename, 'r', 'ieee-le');

status = fseek(fid, offset, 'bof');
if status==-1
  fclose(fid);
  error(message('MATLAB:imagesci:readbmpdata:bad32bitDataOffset'));
end

X = fread(fid, numSamples, 'uint8=>uint8');
fclose(fid);

% Fill in any missing values with zeros.
if (numel(X) < numSamples)
    warning(message('MATLAB:imagesci:readbmpdata:truncatedImageData'));
    X(numSamples) = 0;
end

if (height >= 0)
  X = rot90(reshape(X, scanlineLength, abs(height)));
else
  X = reshape(X, scanlineLength, abs(height))';
end

RGB(1:abs(height), 1:width, 3) = X(:,1:4:end);
RGB(:, :, 2) = X(:,2:4:end);
RGB(:, :, 1) = X(:,3:4:end);
