function [b , preD] = pcaRegression(DIN, SCRS, mth)
%% pcaRegression:
%
%
% Usage:
%   [b , preD] = pcaRegression(DIN, SCRS, mth)
%
% Input:
%   DIN: raw data
%   SCRS: PC Scores
%   mth: method to obtain covariance [pcr|cca]
%
% Output:
%   b:
%   preD: predicted values given input data
%

%% Determine  method for regressino
if nargin < 3
    % Default to PCR
    mth = 'pcr';
end

switch mth
    case 'pcr'
        %% 
        ND   = DIN - mean(DIN);
        b    = SCRS \ ND;
        preD = SCRS * b;
        
    case 'cca'
        %% 
        X           = SCRS;
        Y           = DIN - mean(DIN);
        [A,B,R,U,V] = canoncorr(X,Y);
        
        b    = A;
        preD = U / B;
        
    otherwise
        fprintf(2, 'Regression Method %s must be [pcr|cca]\n', mth);
        [b , preD] = deal([]);
end
end


