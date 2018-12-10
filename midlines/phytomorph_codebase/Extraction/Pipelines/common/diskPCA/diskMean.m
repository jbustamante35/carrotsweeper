%function [U] = diskMean(dataNozzel,toVecFunc,selVecFunc)
function [U] = (dataNozzel)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for each vector source
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = zeros(dataNozzel.size());
    while dataNozzel.hasNext()
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the next vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = dataNozzel.next();
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply func:Tensor-->vector(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = toVecFunc(I,e);
        %}
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply func:vector(s)-->vector(s)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = selVecFunc(I,e);
        %}
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init sum
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if e == 1
            s = zeros(size(I,1),1);
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % accumulate vectors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        s = s + sum(I,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get size for expected value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z(e) = size(I,2);
    end
    % get expected value
    U = s * sum(z)^-1;
end