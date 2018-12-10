function [bTs] = ipca(t)
    bTs = myE();    
    %%%%%%%%%%%
    % get information
    n = t.rank();
    rnk = (n-1);
    %%%%%%%%%%%
    % loop over rank
    for i = 1:rnk        
        % decompose
        [u e info(i).l info(i).pe] = decomp(t,i);        
        % add basis to set
        bTs.addBasis(e,u);        
    end
    bTs.attachInformation(info);
end

function [u e l pe] = decomp(d,dim)
    % cycle permute
    d.cp(dim-1);
    % reshape the tensor
    d.r();        
    % decompose
    [u e l pe] = p(d.d);
    % inverse reshape        
    d.i();  
end

%%%%%%%%%%%%%%%%
% eigen decomposition of d for k
function [u v l pe] = p(d,k)    
    warning off;
    if nargin == 1
        k = size(d,1);
    end
    u = 0;
    %u = mean(d,2);
    %d = bsxfun(@minus,d,u);
    n = size(d,2);
    c = (1/n)*(d*d');
    [v l] = eigs(c,k);
    l = diag(l);
    % compute percent explained
    pe = l*sum(l)^-1;
    pe = cumsum(pe);
end

%%%%%%%%%%%%%%%%
% error measure
function [e] = em(d,s)
    e = sum((d-s).*(d-s),2).^.5;
end

    %{
    %%%%%%%%%%%
    % select basis on threshold of percent explained
    numB = selectBasisNumber(o,.98);    
    o = selectBasis(o,numB);
    % rotate trials to first dim
    t = rotT(d);
    % build out the basis vectors
    bf = buildBF(o);
    bf.d = hilbertNormalize(bf.d')';
    
    % note: NOT SURE ABOUT orthgonal basis
    pck.comp = bf.d'*t;
    % simulate
    sim = bf.d*pck.comp;
    % reshape simulation
    % vsim = reshape(sim,size(d));
    % data and sim error
    pck.error = em(t',sim');
    % construct package
    pck.o = o;
    pck.bf = bf.d;
    pck.comp = pck.comp';
    %}