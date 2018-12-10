function comprout=writegif(mat,map,filename,varargin)
%WRITEGIF Write a GIF (Graphics Interchange Format) file to disk.
%	WRITEGIF(X,MAP,FILENAME) writes the indexed image X,MAP
%   to the file specified by the string FILENAME. The extension 
%   '.gif' will be added to FILENAME if it doesn't already 
%   have an extension.
%
%   X can be a single image (M-by-N) or a series of
%   frames(M-by-N-by-1-by-P).  X must be of type uint8, logical, or double.
%   If X is uint8 or logical, then colormap indexing starts at zero.  If X
%   is double, then colormap indexing starts at one.
%
%   MAP can be a single colormap (M-by-3) that is applied to all frames, or a series of colormaps
%   (M-by-3-by-P) where P is the number of frames specified in X.  
%
%   WRITEGIF(X,[],FILENAME) writes the grayscale image GRAY
%   to the file.
%
%   WRITEGIF(...,'writemode','append') appends a single image to a file existing
%   on disk.
%
%   WRITEGIF(...,'comment',TEXT) writes an image to a file with the
%   comment specified in TEXT.  The comment is placed immediately before the image.
%   TEXT can be either a character array or a cell array of strings.  If
%   TEXT is a cell array of strings, then a carriage return is added after
%   each row.
%
%   WRITEGIF(...,'disposalmethod',DMETH) specifies the disposal
%   method for the image.  DMETH must be one of
%   'leaveInPlace','restoreBG','restorePrevious', or 'doNotSpecify'.
%
%   WRITEGIF(...,'delaytime',TIME) specifies the delay before displaying
%   the next image.  TIME must be a scalar value measured in seconds
%   between 0 and 655 inclusive.
%   
%   WRITEGIF(...,'transparentcolor',COLOR) specifies the transparent color
%   for the image.  COLOR must be a scalar index into the colormap.  If X
%   is uint8 or logical, then indexing starts at 0.  If X is double, then
%   indexing starts at 1.
%
%   WRITEGIF(...,'backgroundcolor',COLOR) specifies the background color
%   for the image.  COLOR must be a scalar index into the colormap.  If X
%   is uint8 or logical, then indexing starts at 0.  If X is double, then
%   indexing starts at 1.
%
%   WRITEGIF(...,'loopcount',COUNT) specifies the number of times to repeat
%   the animation.  If COUNT is Inf, the animation will be continuously
%   looped.  If COUNT is 0, the animation will be played once.  If
%   COUNT is 1, the animation will be played twice, etc.  The maximum value
%   of COUNT (other than Inf) is 65535.

%   WRITEGIF(...,'screensize',SIZE) specifies the screensize for the
%   image.  SIZE must be a two element vector where the first element is
%   the screen height and the second element is the screen width.
%
%   WRITEGIF(...,'location',LOC) specifies the offset of the top left
%   corner of the image relative to the top left corner of the screen.  LOC
%   must be a two element vector where the first element is the offset from
%   the top and the second element is the offset from the left.   
%
%	See also: GIFREAD, BMPWRITE, HDFWRITE, PCXWRITE, TIFFWRITE,
%	          XWDWRITE.

%	Copyright 1993-2010 The MathWorks, Inc.
%	$Revision: 1.1.6.11 $  $Date: 2011/05/17 02:28:13 $

propStrings = {'writemode','comment','disposalmethod',...
        'delaytime','transparentcolor',...
        'loopcount','location','screensize','backgroundcolor'};

%minimum of two arguments
error(nargchk(2, Inf, nargin));

%check filename string
if ~ischar(filename) || isempty(filename)
  error(message('MATLAB:imagesci:writegif:badFilename'));
end

%check for number of dimensions greater than four
nd = ndims(mat);

if nd > 4
  error(message('MATLAB:imagesci:writegif:tooManyDimensions', nd));
end

%RGB output not supported
if (nd == 3)
    if (size(mat,3) == 3)
        error(message('MATLAB:imagesci:writegif:rgbUnsupported'));        
    else
        error(message('MATLAB:imagesci:writegif:tooManyColorPlanes', nd));
    end
end

%check image data type
if (~isa(mat, 'logical') && ~isa(mat,'uint8') && ~isa(mat,'double'))
    error(message('MATLAB:imagesci:writegif:badType'));
end

ndcmap = ndims(map);
if ndcmap > 3
  error(message('MATLAB:imagesci:writegif:badColormapDimensions', ndcmap));
elseif (ndcmap == 3) && ~(nd == 4)
  error(message('MATLAB:imagesci:writegif:badColormapForSingleFrame'));
end

if size(map,1)>256
    warning(message('MATLAB:imagesci:writegif:tooManyColormapEntries'));
    map = map(1:256,:,:);
end
    
%this variable determines whether or not color index params are zero based
%or one based (this is determined by the data type of mat)
zerobased = 1;

%convert to zero based indexing if data is one based (type double)
if(isa(mat,'double'))
    mat = mat-1;
    zerobased = 0;
end

%%%%%%%%%%%%%
%Commented code dealing with userinput and aspectratio
%are commented because the user does not need control
%over these GIF features.
%%%%%%%%%%%%%

%parameter defaults
writemode = 0;  %default writemode is replace (0)
userinput = 0;
disposalmethod = 0;
delaytime = 50;
transparentcolor = -1;
comment = '';  
%bitdepth = ceil(log2(size(map,1)));
backgroundcolor = 0;
aspectratio = 0;
loopcount = -1;  %default value from user's perspective is actually 0
location = [0 0];
screensize = [size(mat,2) size(mat,1)];


% Process param/value pairs
propStrings = {'writemode','comment','disposalmethod',...
        'delaytime','transparentcolor',...
        'loopcount','location','screensize','backgroundcolor'};
%        'userinput','aspectratio'...

%Keep track of which file scope parameters
%are requested so that a warning can be generated
%if writemode is append.
commentrequested = false;
loopcountrequested = false;
screensizerequested = false;
backgroundcolorrequested = false;

for k = 1:2:length(varargin)
    prop = lower(varargin{k});
    value = varargin{k+1};
    if (~ischar(prop))
        error(message('MATLAB:imagesci:writegif:badParameterName'));
    end
    idx = find(strncmp(prop, propStrings, numel(prop)));
    if (isempty(idx))
        error(message('MATLAB:imagesci:writegif:unrecognizedParameter', prop));
    elseif (length(idx) > 1)
        error(message('MATLAB:imagesci:writegif:ambiguousParameter', prop));
    end
    
    prop = deblank(propStrings{idx});
    
    switch prop
        case 'writemode'
            value = lower(value);
            switch value;
                case 'append'
                    writemode = 1;
                case 'overwrite'
                    writemode = 0;
                otherwise
                    error(message('MATLAB:imagesci:writegif:badWriteMode'));
            end
        case 'comment'
            comment = value;
            if (~ischar(comment) & ~iscellstr(comment))
                error(message('MATLAB:imagesci:writegif:badComment'));
            elseif(iscellstr(comment))
                tmp = '';
                for str = comment
                    %TODO - make sure newline characters are correct
                    tmp = [tmp str{1} char(13) char(10)];
                end
                comment = tmp;
            end
            commentrequested = true;
        case 'disposalmethod'
            value = lower(value);
            switch value;
                case 'leaveinplace'
                    disposalmethod=1;
                case 'restorebg'
                    disposalmethod=2;    
                case 'restoreprevious'
                    disposalmethod=3;   
                case 'donotspecify'
                    disposalmethod=0;    
                otherwise
                    error(message('MATLAB:imagesci:writegif:badDisposalMethod'));
            end
        case 'delaytime'
            if(isnumeric(value) && (value<=655 && value>=0))
                %convert seconds to the nearest hundredth of a second
                delaytime = round(value*100);
            else
                error(message('MATLAB:imagesci:writegif:badDelayTime'));
            end
        case 'transparentcolor'
            if(isnumeric(value))
                transparentcolor = round(double(value));
                transparentcolor = transparentcolor - double(not(zerobased));  %if one-based image data, then this is also treated as one-based
            else
                error(message('MATLAB:imagesci:writegif:badTransparentColor'));
            end
        case 'loopcount'
            if(isnumeric(value) && (value==inf || (value<=65535 && value>=0)))
                loopcount = round(value);
                switch(loopcount)
                    case 0
                        loopcount = -1;
                    case Inf
                        loopcount = 0;
                end
            else
                error(message('MATLAB:imagesci:writegif:badLoopCounter'));
            end
            loopcountrequested = true;
        case 'location'
            if(isnumeric(value) && ndims(value)==2 && length(value)==2 && min(size(value))==1)
                location = round(value);
            else
                error(message('MATLAB:imagesci:writegif:badLocation'));
            end
        case 'screensize'
            if(isnumeric(value) && ndims(value)==2 && length(value)==2 && min(size(value))==1)
                screensize = round(value);
            else
                error(message('MATLAB:imagesci:writegif:badScreenSize'));
            end
            screensizerequested=true;
        case 'backgroundcolor'
            if(isnumeric(value))
                backgroundcolor = double(value);
                if(not(zerobased))
                    backgroundcolor = backgroundcolor - 1;  %if one-based image data, then this is also treated as one-based
                end
            end
            backgroundcolorrequested = true;
%         case 'aspectratio'
%             tmp = (value*64)/15;
%             if(isnumeric(value) && value>=0 && value<=255)
%                 aspectratio = round(tmp);
%             else
%                 error('AspectRatio must be an integer between 0 and the 255 inclusive');
%             end
%             if(writemode)
%                 warning('The AspectRatio parameter has no effect when appending to a file');
%             end
%         case 'userinput'
%             switch value;
%                 case 0
%                     userinput=0;
%                 case 1
%                     userinput=1;    
%                 otherwise
%                     error('UserInput must be either 0 or 1');
%             end            
    end
end

%add .gif extension to filename if necessary
if isempty(strfind(filename,'.'))
    filename = [filename '.gif']; 
end

%warn if file scope parameters requested in append mode
if(writemode == 1)
    if(commentrequested)
        warning(message('MATLAB:imagesci:writegif:commentInAppendMode'))
        comment = '';
    end
    if(loopcountrequested)
        warning(message('MATLAB:imagesci:writegif:loopcountInAppendMode'))
        loopcount = -1;
    end
    if(screensizerequested)
        warning(message('MATLAB:imagesci:writegif:screensizeInAppendMode'))
        screensize = [size(mat,2) size(mat,1)];
    end
    if(backgroundcolorrequested)
        warning(message('MATLAB:imagesci:writegif:backgroundcolorInAppendMode'))
        backgroundcolor = 0;
    end
end

%generate grayscale colormap if map is empty and we are not appending 
%(will use global colormap if appending)
if(isempty(map) && not(writemode==1))
    if(islogical(mat))
        map = gray(2);
    else
        map = gray(256);
    end
end

% assume normal orientation of the colormap and image data passed in
if nd==2
    mat = mat';
elseif nd ==4
    mat = permute(mat,[2,1,3,4]);
end

if ndcmap==2
    map = map';
elseif ndcmap==3
    map = permute(map,[2,1,3]);
end

if(isempty(map))
    cmaplength = 256;
else
    cmaplength = size(map,2);
end

greaterthancmap = find(mat>=cmaplength);
lessthancmap = find(mat<0);

if(~isempty(greaterthancmap) || ~isempty(find(mat<0)))
    warning(message('MATLAB:imagesci:writegif:dataOutOfRange'));
    mat(greaterthancmap) = cmaplength-1;
    mat(lessthancmap) = 0;
end

if(transparentcolor>=cmaplength || transparentcolor < -1)
    warning(message('MATLAB:imagesci:writegif:transparentcolorNotUsed'));
    transparentcolor = -1;
end
if(backgroundcolor>=cmaplength || backgroundcolor < -1)
    warning(message('MATLAB:imagesci:writegif:backgroundcolorNotUsed'));
    backgroundcolor = -1;
end

%make sure data is uint8 format
mat = uint8(mat);

wgifc(mat, map, filename,writemode,userinput,disposalmethod,delaytime,...
    transparentcolor,comment,backgroundcolor,aspectratio,loopcount,location,screensize);
