function [ truncSpline, length ] = calcAdjSpikeLength( posAdj, spike  )
%CALCADJSPIKELENGTH Calculates length from an adjusted position.
%   

    %posAdj = find(spike(:,1) == pos);
    trunc = 1:posAdj;
    if ~isempty(posAdj) 
        if posAdj ~= 1
            truncSpline = csaps(trunc, spike(trunc,:)', 10e-8);
        else
            truncSpline = csaps(1:size(spike, 1), spike(:,:)', 10e-8);
        end
    else
        truncSpline = csaps(1:size(spike, 1), spike(:,:)', 10e-8);
    end

    length = calcArcLength(truncSpline);
end

