function [cm,gn] = confusionmat(g,ghat,varargin)
% CONFUSIONMAT Confusion matrix for classification algorithms.
%    CM = CONFUSIONMAT(G,GHAT) returns the confusion matrix CM determined
%    by the known group labels G and the predicted group labels GHAT. G and
%    GHAT are grouping variables with the same number of observations. G
%    and GHAT can be categorical, numeric, or logical vectors;
%    single-column cell arrays of strings; or character matrices (each row
%    representing a group label). G and GHAT must be of the same type. CM
%    is a square matrix with size equal to the total number of distinct
%    elements in G and GHAT. CM(I,J) represents the count of instances
%    whose known group labels are group I and whose predicted group labels
%    are group J. CONFUSIONMAT treats NaNs, empty strings or 'undefined'
%    values in G or GHAT as missing values, and the corresponding
%    observations are not counted.
%
%    The sets of groups and the orders of group labels in rows and
%    columns of CM are the same. They include all the groups appearing in
%    GN, and have the same order of group labels as GN, where GN is the
%    second of output of grp2idx([G;GHAT]).
%
%    CM = CONFUSIONMAT(G,GHAT,'ORDER',ORDER) returns the confusion matrix
%    with the order of rows (and columns) specified by ORDER.  ORDER is a
%    vector containing group labels and whose values can be compared to
%    those in G or GHAT using the equality operator. ORDER must contain all
%    the labels appearing in G or GHAT. ORDER can contain labels which do
%    not appear in G and GHAT, and hence CM will have zeros in the
%    corresponding rows and columns.
%
%    [CM, GORDER] = CONFUSIONMAT(G, GHAT) returns the order of group labels
%    for rows and columns of CM. GORDER has the same type as G and GHAT.
%
%   Example:
%      % Compute the resubstitution confusion matrix for applying CLASSIFY
%      % on Fisher iris data.
%      load fisheriris
%      x = meas;
%      y = species;
%      yhat = classify(x,x,y);
%      [cm,order] = confusionmat(y,yhat);
%
%   See also CROSSTAB, GRP2IDX.

%   Copyright 2008-2016 The MathWorks, Inc.


if nargin < 2
    error(message('stats:confusionmat:NotEnoughInputs'));
end

gClass = class(g);
if ~strcmp(gClass,class(ghat))
    error(message('stats:confusionmat:GTypeMismatch'));
end

if ~isnumeric(g) && ~islogical(g) && ~isa(g,'categorical') ...
        && ~iscellstr(g)  && ~ischar(g) && ~isstring(g)
    error(message('stats:confusionmat:GTypeIncorrect'));
end

if ischar(g)
    if ~ismatrix(g) || ~ismatrix(ghat)
        error(message('stats:confusionmat:BadGroup'));
    end
    g = cellstr(g);
    ghat = cellstr(ghat);
elseif ~isvector(g) || ~isvector(ghat) 
    error(message('stats:confusionmat:BadGroup'));
else
    g = g(:);
    ghat = ghat(:);
end


if size(g,1) ~= size(ghat,1)
    error(message('stats:confusionmat:GRowNumMismatch'));
end

if isa(g,'categorical')
    if isordinal(g)
        if ~isordinal(ghat)
            error(message('stats:confusionmat:GOrdinalLevelsMismatch'));
        elseif ~isequal(categories(g),categories(ghat))
            error(message('stats:confusionmat:GOrdinalLevelsMismatch'));
        end
    elseif isordinal(ghat)
        error(message('stats:confusionmat:GOrdinalLevelsMismatch'));
    end
end

pnames = {'order'};
dflts =  {[]};
[order] = internal.stats.parseArgs(pnames, dflts, varargin{:});

if ~isempty(order)
    if ischar(order)
        if ~ismatrix(order)
            error(message('stats:confusionmat:NDCharArrayORDER'));
        end
        order = cellstr(order);
    elseif ~isvector(order)
        error(message('stats:confusionmat:NonVectorORDER'));
    end

    if isa(g,'categorical')
        if iscellstr(order)
            if any(strcmp('',strtrim(order)))
                error(message('stats:confusionmat:OrderHasEmptyString'));
            end
        elseif iscategorical(order)
            if any(isundefined(order))
                error(message('stats:confusionmat:OrderHasUndefined'));
            end
        else
            error(message('stats:confusionmat:TypeMismatchOrder'));
        end
        
        g = setcatsLocal(g,order);
        ghat = setcatsLocal(ghat,order);
        
    else % g is not categorical vector

        if isnumeric(g)
            if ~isnumeric(order)
                error(message('stats:confusionmat:TypeMismatchOrder'));
            end
            if any(isnan(order))
                error(message('stats:confusionmat:OrderHasNaN'));
            end

        elseif islogical(g)
            if islogical(order)
                %OK. do nothing
            elseif isnumeric(order)
                if any(isnan(order))
                    error(message('stats:confusionmat:OrderHasNaN'));
                end

                order = logical(order);
                
            else
                error(message('stats:confusionmat:TypeMismatchOrder'));
            end

        elseif iscellstr(g)
            if ~iscellstr(order)
                error(message('stats:confusionmat:TypeMismatchLevels'));
            end
            if any(strcmp('',strtrim(order)))
                error(message('stats:confusionmat:OrderHasEmptyString'));
            end
        end

        try
            uorder = unique(order);
        catch ME
            error(message('stats:confusionmat:UniqueMethodFailedOrder'));
        end

        if length(uorder) < length(order)
            error(message('stats:confusionmat:DuplicatedOrder'));
        end

        order = order(:);
    end
end

gLen = size(g,1);
[idx,gn,gl] = grp2idx([g;ghat]);
gidx = idx(1:gLen);
ghatidx = idx(gLen+1:gLen*2);

%ignore NaN values in GIDX and GHATIDX
nanrows = isnan(gidx) | isnan(ghatidx);
if any(nanrows)
    gidx(nanrows) = [];
    ghatidx(nanrows) = [];
end

gLen = size(gidx,1);
gnLen = length(gn);

cm = accumarray([gidx,ghatidx], ones(gLen,1),[gnLen, gnLen]);

%convert gn to the same type as g
gn = gl;

if ~isempty(order)
    if ~isa(g,'categorical')
        %get the map from the default order to the given order
        [hasAllLabel,map] = ismember(gn,order);
        
        if ~all(hasAllLabel)
            error(message('stats:confusionmat:OrderInsufficientLabels'));
        end
        orderLen = length(order);
        cm2 = zeros(orderLen, orderLen);
        cm2(map,map) = cm(:,:);
        cm = cm2;

        if nargout > 1
            %convert gn to the same type as g
            if strcmp(gClass,'char')
                gn = char(order);
            else
                gn = order;
            end
        end
    end
elseif strcmp(gClass,'char')
    gn = char(gn);
end


function b = setcatsLocal(a,newCategories)
if iscategorical(newCategories)
    if ~isa(newCategories,class(a))
        error(message('stats:confusionmat:TypeMismatchOrder'));
    elseif isordinal(a)
        if ~isordinal(newCategories)
            error(message('stats:confusionmat:TypeMismatchOrder'));
        elseif ~isequal(categories(a),categories(newCategories))
            error(message('stats:confusionmat:TypeMismatchOrder'));
        end
    elseif isordinal(newCategories)
        error(message('stats:confusionmat:TypeMismatchOrder'));
    end
    newCategories = cellstr(newCategories);
end

existingCategories = categories(a);
if ~isempty(setdiff(existingCategories,newCategories))
    error(message('stats:confusionmat:OrderInsufficientLabels'));
end

b = addcats(a,newCategories,'after',existingCategories{end});
b = reordercats(b,newCategories);
