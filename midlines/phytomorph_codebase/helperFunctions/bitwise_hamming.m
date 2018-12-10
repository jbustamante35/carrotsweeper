%BITWISE_HAMMING Compute all hamming distances between two sets of bit vectors
%
%   C = bitwise_hamming(A, B)
%
% Given two sets of bit vectors (each column being a bit vector), compute
% the hamming distances between all pairs of vectors between the two sets.
%
%IN:
%   A - MxP matrix of bit vectors.
%   B - MxQ matrix of bit vectors.
%
%OUT:
%   C - PxQ matrix of bitwise hamming distances.

function C = bitwise_hamming(A, B)

% Get the input sizes
[a, a] = size(A);
[b, b] = size(B);

% Typecast to uint8
A = reshape(typecast(A(:), 'uint8'), [], a);
B = reshape(typecast(B(:), 'uint8'), [], b);
c = size(B, 1);
assert(size(A, 1) == c, 'Input columns must be of the same bit-length');

% Check if joint array less than 64MB
if numel(A) * b < 67108864
    % Do the bitxor in one go
    A = repmat(A, [1 1 b]);
    B = repmat(reshape(B, c, 1, b), [1 a 1]);
    C = reshape(bitcount(bitxor(A, B)), a, b);
else
    % Do the bitxor in a loop, to avoid too much memory use
    if b > a % Do the smallest for loop
        for i = a:-1:1
            C(i,:) = bitcount(bitxor(repmat(A(:,i), [1 b]), B));
        end
    else
        for i = b:-1:1
            C(:,i) = bitcount(bitxor(repmat(B(:,i), [1 a]), A))';
        end
    end
end
end
