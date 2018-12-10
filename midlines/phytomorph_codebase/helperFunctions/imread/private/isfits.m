function tf = isfits(filename)
%ISFITS Returns true for a FITS file.
%   TF = ISFITS(FILENAME)

%   Copyright 1984-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2005/11/15 01:12:31 $

fid = fopen(filename, 'r');
if (fid < 0)
    tf = false;
else
    sig = fread(fid, 6, 'char=>char')';
    fclose(fid);
    tf = isequal(sig, 'SIMPLE');
end
