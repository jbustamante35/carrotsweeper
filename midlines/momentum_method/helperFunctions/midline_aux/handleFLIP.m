function [I direc] = handleFLIP(I,direc)
if isempty(direc)
    % handle the direction
    [direc] = DC(I);
end
switch direc
   case 1                                                      % from the left
       % do nothing
   case 3                                                      % from the right
       I = fliplr(I);                   
   case 2                                                      % from the top
       I = I';
   case 4                                                      % from the bottom
       I = fliplr(I');
end
% handle the direction   
% read and flip
% fliped the right and top on november 4 2015