function data = h5postprocessopaques(~,space,raw_data)
% Opaque data is a sequence of unstructured bytes where the only piece
% of information is the number of bytes per element.  The output here will
% be a cell array with each cell element being an nx1 column of uint8 data.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/10/11 14:50:57 $


[ndims,h5_dims] = H5S.get_simple_extent_dims(space);
dims = fliplr(h5_dims);

if ndims == 1
    % The dataset is one-dimensional.  Force the output to be a column.
    data = cell(dims(1),1);
else
    data = cell(dims);
end
for j = 1:numel(data)
    data{j} = raw_data(:,j);
end
