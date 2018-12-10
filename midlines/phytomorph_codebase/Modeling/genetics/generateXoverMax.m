function [xM] = generateXoverMax(sites,len,exeThres,type)
    % generate rand for threshold compare -- better later
    exeProb = rand(1,numel(exeThres));
    toExe = exeProb > exeThres;
    sites = sites(toExe);
    % stack of potential sites
    sites = [1 sites len];
    % fill value
    b = logical(1);
    switch type
        case 'full'
            % init the diag of xvoer matrix
            xM = zeros(1,len);
            % for each haplotype
            for e = 1:(numel(sites)-1)
                % fill in the diag
                xM(sites(e):sites(e+1)) = b;
                % flip for next fill - assume bi-stable state
                b = ~b;
            end
            xM = [[diag(xM) diag(~xM)];[diag(~xM) diag(xM)]];
        case 'sparse'
            x = 1:len;
            y = 1:len;
            % init the diag of xvoer matrix
            xMd = zeros(1,len);
            % for each haplotype
            for e = 1:(numel(sites)-1)
                % fill in the diag
                xMd(sites(e):sites(e+1)) = b;
                % flip for next fill - assume bi-stable state
                b = ~b;
            end
            xM = sparse(x,y,xMd,len,len);
            nxM = sparse(x,y,~xMd,len,len);
            xM = [[xM nxM];[nxM xM]];
end