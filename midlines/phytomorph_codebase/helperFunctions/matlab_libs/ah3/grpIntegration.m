function [v vi] = grpIntegration(gVec,data,target)
    %%%%%%%%%%%%%%%%%%%%%%%%
    % gVec      := the group vector labels
    % data      := the data vectors - gradings on the data
    % target    := the number of expected objects in the group
    %%%%%%%%%%%%%%%%%%%%%%%%
    % find the groups via unique
    UQ = unique(gVec(:,1));
    vi = [];
    %%%%%%%%%%%%%%%%%%%%%%%%
    % loop over all classes
    for u = 1:numel(UQ)
        % find the data from the uth group
        sidx = gVec(:,1) == UQ(u);
        % select the data from the uth group - prob data
        sub_data = data(sidx);
        % if nargout request > 2
        if (nargout==2)
            % get the prob of the top n
            [sig_data sidx] = sort(sub_data);
            si = sidx(1:target);
            z = zeros(size(sidx));
            vi = [vi;z(si)];
        end
        % sum the prob
        v(u) = sum(sub_data);
    end
end