function [data,map] = hdfraster8read(hinfo)
%HDFRASTER8READ
%
%   [DATA,MAP] = HDFRASTER8READ(HINFO) returns in the variable DATA the
%   image from the file for the particular 8-bit raster image described by
%   HINFO.  MAP contains the colormap if one exists for the image.  HINFO is
%   A structure extracted from the output structure of HDFINFO.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/05/17 02:27:34 $

parseInputs(hinfo);

status = hdfdfr8('readref',hinfo.Filename,hinfo.Ref);
if status == -1
    hdferrmsg(status,'DFR8','restart');
end

[data,map,status]  = hdfdfr8('getimage',hinfo.Filename);
if status == -1
    hdferrmsg(status,'DFR8','getimage');
end

status = hdfdfr8('restart');
if status == -1
    hdferrmsg(status,'DFR8','restart');
end

%Put the image data and colormap in the right order for image display in
%MATLAB
data = data';
map = double(map')/255;
return;

%=======================================================================
function parseInputs(hinfo)

error(nargchk(1,1,nargin, 'struct'));

%Verify required fields

if ~isstruct(hinfo)
  error(message('MATLAB:imagesci:hdfraster8read:notStruct'));
end
fNames = fieldnames(hinfo);
numFields = length(fNames);
reqFields = {'Filename','Ref'};
numReqFields = length(reqFields);
if numFields >= numReqFields
  for i=1:numReqFields
    if ~isfield(hinfo,reqFields{i})
      error(message('MATLAB:imagesci:hdfraster8read:missingRequiredField'));
    end
  end
else 
  error(message('MATLAB:imagesci:hdfraster8read:invalidInputs'));
end
return;





