function [sam E_store PL] = bugEye(I,NH,PL)
    if nargin <= 2
        E = edge(I);
        E = imclearborder(E,8);
        [e2 e1] = find(E);
        PL = [e1 e2];
    end
    % gradient of image
    [dx dy] = gradient(I);
    sam = zeros([size(NH,2) size(PL,1)]);
    for e = 1:size(PL,1)
        % at point
        pt = [PL(e,1);PL(e,2)];
        % normal space at the point is along the gradient
        nor = [dy(PL(e,2),PL(e,1));dx(PL(e,2),PL(e,1))];
        % normalize the normal vector
        nor = nor / norm(nor);
        % create the basis space
        E = [nor [-nor(2);nor(1)]];
        % rotate the nhood
        NHt = E'*NH;
        % move the nhood
        NHt = bsxfun(@plus,NHt,pt);
        % interp1 the nhood post rotate and move
        sam(:,e) = ba_interp2(I,NHt(1,:)',NHt(2,:)');
        if nargout >=2
            E_store(:,:,e) = E;
        end
    end
end