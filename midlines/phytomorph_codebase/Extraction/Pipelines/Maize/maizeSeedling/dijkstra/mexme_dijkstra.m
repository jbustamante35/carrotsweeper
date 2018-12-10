try
    C   = computer;
    
    if(isempty(strfind(C , '64')))
        echo on
        mex dijkstra.c
        mex Kadjacency.c
        mex Radjacency.c
        echo off
        
    else
        echo on
        mex dijkstra.c -largeArrayDims
        mex Kadjacency.c -largeArrayDims
        mex Radjacency.c -largeArrayDims
        echo off
    end
catch exception
    if(~isempty(exception))
        fprintf(['\n Error during compilation, be sure to:\n'...
            'i)  You have C compiler installed (prefered compiler are MSVC/Intel/GCC)\n'...
            'ii) You did "mex -setup" in matlab prompt before running mexme_dijkstra\n']);
    end
end

