function [strain] = convertToStrain(sig,filter,tau)
    ssig = imfilter(sig,fspecial('average',[1 filter]),'replicate');
    
    strain = zscore(ssig,1,2);
    
    %{
    init = mean(ssig(:,1:min(size(ssig,2),tau)),2);
    %final = mean(ssig(:,(end-tau):end),2);
    ssig = bsxfun(@minus,ssig,init);
    srt = sort(ssig,2,'descend');
    mx = mean(srt(:,1:2),2);
    mn = mean(srt(:,(end-2):end),2);
    %final = mean(ssig(:,(end-tau):end),2);
    upper = ssig.*(ssig > 0);
    lower = ssig.*(ssig < 0);
    mx(mx < 1) = 1;
    mn(mn > -1) = 1;
    upper = bsxfun(@times,upper,mx.^-1);
    lower = bsxfun(@times,lower,-mn.^-1);
    strain = upper + lower;
    %strain = bsxfun(@times,ssig,mx.^-1);
    %}
end