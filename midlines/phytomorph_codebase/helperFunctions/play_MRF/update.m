function [snp0 snp1] = update(p0,p1,MTI,midx0,midx1,SZ,p_O_0,p_O_1)

    p0(p0==0) = eps;
    p1(p1==0) = eps;
    T = p0 + p1;
    p1 = p1 .* T .^-1;
    p0 = p0 .* T .^-1;

    lp0 = log(p0);
    lp1 = log(p1);

    %lp0 = padarray(lp0, [1 1], 'replicate', 'both');
    %lp1 = padarray(lp1, [1 1], 'replicate', 'both');
    
    rp0 = im2colF(lp0,[3 3],[1 1]);
    rp1 = im2colF(lp1,[3 3],[1 1]);
    newI = mtimesx(double(MTI),[rp1;rp0]);
    newI0 = newI(midx0,:);
    newI1 = newI(midx1,:);
    
    
    %newI0 = mtimesx((P(midx0)+1),'T',newI0);
    %newI1 = mtimesx(P(midx1),'T',newI1);
    %newI0 = mean(newI0,1)*sum(midx0);
    %newI1 = mean(newI1,1)*sum(midx1);
    
    
    newI1 = max(newI1,[],1);
    %newI0 = max(newI0,1);
    
    newI0 = log(1 - exp(newI1));
    R = newI1 - newI0;
    R = (exp(R)+1).^-1;
    
    v0 = reshape(R,SZ);
    v1 = 1-v0;
    %{
    v1 = reshape(exp(A*newI1),SZ);
    v0 = reshape(exp(A*newI0),SZ);
    %}
    T = v1 + v0;
    v1 = v1 .* T.^-1;
    v0 = v0 .* T.^-1;

    np1 = padarray(v1, [1 1], 'replicate', 'both');
    np0 = padarray(v0, [1 1], 'replicate', 'both');
    
    
    snp0 = np0.^1.*p_O_0;
    snp1 = np1.^1.*p_O_1;
    
    
    TN = snp1 + snp0;
    snp1 = snp1 .* TN.^-1;
    snp0 = snp0 .* TN.^-1;
    
    
    snp0(isnan(snp0(:))) = .5;
    snp1(isnan(snp1(:))) = .5;

    

end