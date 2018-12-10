function data = hdfraster24read(hinfo)
%HDFRASTER24READ
%
%   DATA = HDFRASTER24READ(HINFO) returns in the variable DATA the image
%   from the file for the particular 24-bit raster image described by HINFO.
%   HINFO is a structure extracted from the output structure of HDFINFO.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/05/17 02:27:32 $

parseInputs(hinfo);
status = hdfdf24('readref',hinfo.Filename,hinfo.Ref);
if status == -1
    hdferrmsg(status,'DF24','readref');
end
	
[data, status] = hdfdf24('getimage',hinfo.Filename);
if status == -1
    hdferrmsg(status,'DF24','getimage');
end
	
status = hdfdf24('restart');
if status == -1
    hdferrmsg(status,'DF24','restart');
end

%Put the image data in the right order for image display in MATLAB
data = permute(data,[3 2 1]);
return;

%=======================================================================
function parseInputs(hinfo)

error(nargchk(1,1,nargin, 'struct'));
	  
%Verify required fields

if ~isstruct(hinfo)
    error(message('MATLAB:imagesci:hdfraster24read:notStruct'));
end
fNames = fieldnames(hinfo);
numFields = length(fNames);
reqFields = {'Filename','Ref'};
numReqFields = length(reqFields);
if numFields >= numReqFields
  for i=1:numReqFields
    if ~isfield(hinfo,reqFields{i})
      error(message('MATLAB:imagesci:hdfraster24read:structMissingRequiredField'));
    end
  end
else 
  error(message('MATLAB:imagesci:hdfraster24read:invalidInputs'));
end
