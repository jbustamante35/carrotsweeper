function [delta,minJump] = phD(baseTheta,theta2,baseFrequency)
    deltaTheta = (theta2 - baseTheta);
    choiceTheta = -sign(deltaTheta).*2*pi + sign(deltaTheta).*deltaTheta;
    choiceVector = [deltaTheta;choiceTheta];
    [~,idx] = min(abs(choiceVector),[],1);
    
    for e = 1:numel(idx)
        delta(e) = choiceVector(idx(e),e);
    end
    
    minJump = sign(delta).*mod(abs(delta),abs(baseFrequency));
    jumpFraction = minJump .* baseFrequency.^-1
    jumpFractionP = jumpFraction - .5
    jumpChoice = [minJump;baseFrequency-minJump];
    
end

%{
    [delta,minJump] = phD(-pi/4,pi-pi/16,2*pi/3);
    [delta,minJump] = phD(28/(40*pi),(8+28)/(40*pi),[3 5]);
    [delta,minJump] = phD(-pi/4,[pi-pi/16 pi/2],[2*pi/3 3*pi/6]);
%}