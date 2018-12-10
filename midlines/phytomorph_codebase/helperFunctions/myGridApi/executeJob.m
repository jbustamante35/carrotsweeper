function [] = executeJob(jP)
    % loop over all jobs
    for j = 1:numel(jP)
        % call each function in the job package
        for fn = 1:numel(jP(j).func)
            jP(j).func(fn).out = jP(j).func(fn).f(jP(j));
        end
    end
end