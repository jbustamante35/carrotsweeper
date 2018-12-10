function [ tort ] = calcTort( S )
%CALCTORT Calculates tortuosity of S
%   Tortuosity calculated as the euclidean distance between ends of S 
%   divided by the line integral of S to give a number between 0 and 1.
%   S: pp spline

    coords = ppval(S, [1 size(S.coefs, 1)/2 ]);
    eucl = sqrt(sum((coords(:,1) - coords(:,2)).^2));

    tort = eucl / calcArcLength(S);


end

