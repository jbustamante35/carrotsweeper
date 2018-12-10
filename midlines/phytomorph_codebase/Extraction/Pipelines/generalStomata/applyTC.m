function [dE,C,errDis] = applyTC(data,U1,E1,U2,E2,modelStruct)
    
    [C,totE,recon,errDis] = drc(data,U1,E1,U2,E2);

       
    totE
    errDis = errDis(:);
    C = [C(:);totE];




    V = C(:);
    if ~isempty(modelStruct)
        for e = 1:numel(V)
            dE(e) = interp1(modelStruct(e).X,modelStruct(e).F,V(e));
        end
        dE = dE';
        dE(end) = exp(-V(end)/25);
    else
        dE = 0;
    end
end


function [C,totE,recon,errDis] = drc(data,U1,E1,U2,E2)


    dsz = size(data);
    %data = zscore(data(:),1,1);
    data = reshape(data,dsz);
    % sim radial vectors
    [C] = PCA_REPROJ_T(data,E1,U1);




    % sim radial vectors
    C = permute(C,[2 1 3]);
    [C] = PCA_REPROJ_T(C,E2,U2);
    recon = PCA_BKPROJ_T(C,E2,U2);
    C = ipermute(C,[2 1 3]);
    recon = ipermute(recon,[2 1 3]);
    recon = PCA_BKPROJ_T(recon,E1,U1);

    

    errDis = recon(:)-data(:);
    totE = norm(errDis);
    errDis = reshape(errDis,size(data));
end