function [out] = mat2pct(in,classType)
    import phytoG.locked.BdataObjects.geometry.implementations.*;    
    import phytoG.locked.BdataStructures.implementations.directMongo.nan.*;
    import phytoG.locked.BdataStructures.implementations.directMongo.gamma.*;
    import java.util.HashMap;
    import java.lang.Double;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create object
    eval(['out = ' classType '();']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % excite to transistion datanode
    % to embedded state
    dE = HashMap();
    dE.put('_dnode',java.lang.Integer(2));
    out.exciteState(dE);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store size
    sz = size(in);
    psz = B2DBL();
    for e = 1:numel(sz)
        psz.add(sz(e));
    end
    out.setProp('size',psz);
    %out.setSize(psz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % store rank
    out.setProp('rank',java.lang.Integer(ndims(in)));
    %out.setRank(java.lang.Integer(ndims(in)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set melted size
    out.setProp('sizeMelted',java.lang.Integer(prod(sz)));
    %out.setSizeM(java.lang.Integer(prod(sz)));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set data
    oVec = javaArray('java.lang.Double',prod(sz));
    for e = 1:numel(in)
        oVec(e) = java.lang.Double(in(e));
    end
    out.setData(oVec);    
end