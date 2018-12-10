function [featureTensor] = fftPatch(I,positionIndex,samplerFunction,freqToKeep)
    patch = samplerFunction(I,positionIndex);
    patch = bsxfun(@minus,patch,mean(patch,2));
    fftResult = fft(patch,[],2);
    ang = unwrap(angle(fftResult),[],2);
    tmpA = ang;
    tmpA = diff(cat(2,tmpA,tmpA(:,end)),1,2);
    tmpA = abs(tmpA)/(2*pi);
    tmpA = .5*(((tmpA.^2)+1).^.5 + (((1-tmpA).^2)+1).^.5);
    patch = abs(fftResult);
    patch = patch(:,freqToKeep);
    tmpA = tmpA(:,freqToKeep);
    featureTensor = cat(3,patch,tmpA);
    
end