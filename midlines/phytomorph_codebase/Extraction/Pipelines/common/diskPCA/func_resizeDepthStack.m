function [I] = func_resizeDepthStack(I,sz,toT,toSort)
    % function to load into nozzle. this function will make thumbnail and
    % sort if toSort flag is true and will stack the color panes
    % horizontally
    % I: image to operate on
    % sz: thumbnail size
    % toT: flag for transpose after thumbnail and sort
    % % create
    % make thumbnail
    I = func_thumbNail(I,sz,false,false);
    % sort if neede
    if toSort
        I = sort(I,1);
    end
    % depth stack the images
    I = func_depthStack(I,toT);
end