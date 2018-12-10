%BITCOUNT Count the number of set bits in each column of the input
%
%   B = bitcount(A)
%
% Count the number of set bits in each column of the input array,
% typecast as a bit vector.
%
%IN:
%   A - MxNx... input array.
%
%OUT:
%   B - 1xNx... output array of bit counts.

function A = bitcount(A)

persistent lutable
if isempty(lutable)
    % Generate the lookup table
    lutable = uint8(sum(dec2bin(0:255) - '0', 2));
end

% Convert to an index into the lookup table
sz = size(A);
sz(1) = 1;
A = reshape(typecast(A(:), 'uint8'), [], prod(sz));

% Look up the number of set bits for each byte
try
    A = intlut(A, lutable);
catch
    A = uint16(A) + uint16(1);
    A = lutable(A);
end

% Sum the number of set bits per column
if size(A, 1) < 32
    A = sum(A, 1, 'native');
else
    A = sum(A, 1);
end
A = reshape(A, sz);
end
