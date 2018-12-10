function [G] = generatePaths(G,maxN,lengthV,N)
    for e = 1:numel(lengthV)
        pl = lengthV(e);
        % create a random number of paths of length pl
        numPaths = randi(maxN,1);
        for n = 1:numPaths
            % select start point from N points
            p1 = randi(N,1);

            % make a chain of length pl
            for e = 1:pl
                % select send point
                p2 = randi(N,1);
                ed = phE([p1 p2]);
                G.addEdge(ed);
                % make end point the next start point
                p1 = p2;
            end
        end
    end
end