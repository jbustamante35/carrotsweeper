function writejp2k(data, map, filename, fmt, varargin)
%WRITEJP2K Internal function facilitating JPEG2000 writes for J2C and JP2.

%   Copyright 2009-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2011/05/17 02:28:15 $

% Input checking.
if (ndims(data) > 3)
    error(message('MATLAB:imagesci:writejp2k:tooManyDims', ndims( data )));
end

if (~isempty(map))
    error(message('MATLAB:imagesci:writejp2k:tooManyDimsForIndexed', ndims( data )));
end

if isfloat(data)
    % single/double data is converted to uint8
    maxval = 255;
    data = uint8(maxval * data);
end

props = set_jp2c_props(data,fmt,varargin{:});

writejp2c(data, filename, props);

function props = set_jp2c_props(data,fmt,varargin)
% SET_JP2C_PROPS
%
% Parse input parameters to produce a properties structure.  
%

%
% Set the default properties.
props.cratio = 1;
props.mode = 'lossy';
props.porder = 'lrcp';
props.qlayers = 1;
props.rlevels = -1;
props.tilewidth = size(data, 2);
props.tileheight = size(data, 1);
props.comment = {};
props.format = fmt;

% Process param/value pairs
paramStrings = {'format'
                'compressionratio'
                'mode'
                'progressionorder'
                'qualitylayers'
                'reductionlevels'
                'tilesize'
                'comment'};

for k = 1:2:length(varargin)
  
    param = lower(varargin{k});
    if (~ischar(param))
        error(message('MATLAB:imagesci:writejp2k:badParameterName'));
    end
    
    idx = find(strncmp(param, paramStrings, numel(param)));
    if (isempty(idx))
        error(message('MATLAB:imagesci:writejp2k:unrecognizedParameter', param));
    elseif (length(idx) > 1)
        error(message('MATLAB:imagesci:writejp2k:ambiguousParameter', param));
    end

    param = deblank(paramStrings{idx});

    props = process_argument_value ( props, param, varargin{k+1} );
    
end

return



% Process a parameter name/value pair, return the new property structure
function output_props = process_argument_value ( props, param_name, param_value )

output_props = props;

switch param_name 
case 'compressionratio'
  
    cratio = param_value;
    
    if (~isa(cratio,'numeric'))
        error(message('MATLAB:imagesci:writejp2k:nonNumericCompressionRatio'))
    end
    if ((cratio < 1) || ~isfinite(cratio))
        error(message('MATLAB:imagesci:writejp2k:badCompressionRatio'));
    end
    
    output_props.cratio = cratio;
    
case 'mode'
    
    mode = lower(param_value);
    
    if ((~ischar(mode)) || ...
        ((~isequal(mode, 'lossy')) && (~isequal(mode, 'lossless'))))
        error(message('MATLAB:imagesci:writejp2k:badMode'))
    end
    
    output_props.mode = mode;
 
case 'progressionorder'
    porder = lower(param_value);
    
    if ((~ischar(porder)) || ...
        (~isequal(porder, 'lrcp') && ~isequal(porder, 'rlcp') && ...
         ~isequal(porder, 'rpcl') && ~isequal(porder, 'pcrl') && ...
         ~isequal(porder, 'cprl')))
        error(message('MATLAB:imagesci:writejp2k:badProgressionOrder'))
    end
    
    output_props.porder = porder;
  
case 'qualitylayers'
    qlayers = param_value;
    
    if (~isa(qlayers,'numeric'))
        error(message('MATLAB:imagesci:writejp2k:nonNumericQualityLayers'))
    end
    if ((qlayers < 1) || (qlayers > 20) || rem(qlayers, 1) ~= 0)
        error(message('MATLAB:imagesci:writejp2k:badQualityLayers'));
    end
    
    output_props.qlayers = qlayers;
    
case 'reductionlevels'
    rlevels = param_value;
    
    if (~isa(rlevels,'numeric'))
        error(message('MATLAB:imagesci:writejp2k:nonNumericReductionLevels'))
    end
    if ((rlevels < 1) || (rlevels > 8)  || rem(rlevels, 1) ~= 0)
        error(message('MATLAB:imagesci:writejp2k:badReductionLevels'));
    end
    
    output_props.rlevels = rlevels;
    
case 'tilesize'
    tilesize = param_value;
    
    if ~isa(tilesize,'numeric') || numel(tilesize) ~= 2
        error(message('MATLAB:imagesci:writejp2k:tileSizeNotTwoElementNumeric'))
    end
    if tilesize(1) < 128 || tilesize(2) < 128 || ...
       ~all(tilesize <= intmax) || any(rem(tilesize, 1) ~= 0)
        error(message('MATLAB:imagesci:writejp2k:badTileSize'));
    end

    output_props.tileheight = tilesize(1);
    output_props.tilewidth = tilesize(2);

case 'comment'
    comment = param_value;
    if (~ischar(comment) && ~iscellstr(comment))
        error(message('MATLAB:imagesci:writejp2k:badComment'));
    end
    % Convert the char matrix to a cell array
    output_props.comment = cellstr(comment);
end

return


