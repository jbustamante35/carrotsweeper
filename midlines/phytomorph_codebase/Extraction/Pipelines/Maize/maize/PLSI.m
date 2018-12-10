function [t p u q w b] = PLSI (x,y)
% 
% Inputs:
% x     x matrix
% y     y matrix
% 
% Outputs:
% t     score for x
% p     loading for x
% u     score for y
% q     loading for y
% b     regression coefficient

% calculate the size of x and y
[nX,mX]  =  size(x);
[nY,mY]  =  size(y);

nMaxIteration = max([mX,mY]);
nMaxOuter = 10000;

for iIteration = 1 : nMaxIteration
    % choose the column of x has the largest square of sum as t.
    % choose the column of y has the largest square of sum as u.    
    [dummy,tNum] =  max(diag(x'*x));
    [dummy,uNum] =  max(diag(y'*y));
    tTemp = x(:,tNum);
    uTemp = y(:,uNum);
    % iteration for outer modeling
    for iOuter = 1 : nMaxOuter
        wTemp = x' * uTemp/ norm(x' * uTemp);
        tNew = x * wTemp;
        qTemp = y' * tNew/ norm (y' * tNew);
        uTemp = y * qTemp;
        if norm(tTemp - tNew) < 10e-15
            break
        end
       
        tTemp = tNew;
        
    end
    % residual deflation:
    bTemp = uTemp'*tTemp/(tTemp'*tTemp);
    pTemp = x' * tTemp/(tTemp' * tTemp);
    x = x - tTemp * pTemp';
    y = y - bTemp * tTemp * qTemp';
    
    % save iteration results to outputs:
    t(:, iIteration)           = tTemp;
    p(:, iIteration)           = pTemp;
    u(:, iIteration)           = uTemp;
    q(:, iIteration)           = qTemp;
    w(:, iIteration)           = wTemp;
    b(iIteration,iIteration)   = bTemp;
    % check for residual to see if we want to continue:
    if (norm(x) ==0 | norm(y) ==0)
        break
    end
end
