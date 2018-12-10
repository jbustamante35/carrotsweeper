function [ret] = site_process_caller(image,site_generate_function,site_process_function,site_reduce_function)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % site_generate_function.func = function to generate a set of sites
    % site_generate_function.OP = options passed to the site generator
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % site_process_function.block_size = size of the blocks to process
    % site_process_function.func = function to process the blocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    % call the function which will generate a set of sites
    fprintf(['Extracting sites - START \n']);tic;
    [sites auxData] = site_generate_function.func(image,site_generate_function.OP);
    fprintf(['Extracting sites - ' num2str(toc) ' -END \n']);
    % loop over the sites and call the block process function with the
    % block size and the image
    for e = 1:numel(sites)        
        ret{e}.data = siteProcess(sites{e},site_process_function.block_size,site_process_function.func,image);
        ret{e}.data = site_reduce_function(ret{e}.data);
        ret{e}.sites = sites{e};
        ret{e}.auxData = auxData(e);
    end    
end