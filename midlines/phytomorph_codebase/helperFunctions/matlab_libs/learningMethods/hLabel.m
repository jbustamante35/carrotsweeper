function [outLabel] = hLabel(vec,inLabel,func)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % notes: added(13.01.25) - this will classify 
    %        but NOT yet determine if the classification
    %        is stable OR "good"   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT:
    %       vec:        vector to recursively group - [trials vector_dim]
    %       inLabel:    init vector labels
    %       func:       set of function handles for labeling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %       outLabel:   output labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % init the out labels
    outLabel = inLabel;    
    % recursive call to label funcs
    for n = 1:numel(func)
        outLabel = f0(vec,outLabel,func{n}.phi,func{n}.para);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iterate over ith group and split into Ni groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outLabel] = f0(vec,inLabel,func,para)
    %%%%%%%%%%%%%%%%%
    % init curLabel to zero
    curLabel = 0;
    outLabel = zeros(size(inLabel));
    UQ0 = unique(inLabel);
    %%%%%%%%%%%%%%%%%
    % for each unique group
    for i = 1:numel(UQ0)
        
        %%%%%%%%%%%%%%%%%
        % isolate the ith group
        idx0 = find(inLabel==UQ0(i));
        %%%%%%%%%%%%%%%%%
        % apply function to the vectors
        %Is = func.filter(vec);
        %%%%%%%%%%%%%%%%%
        % get the ith sub-sample 
        vec_sub = vec(idx0,:);
        %%%%%%%%%%%%%%%%%
        % call to divide the ith sample into groups
        labelVec = func(vec_sub,para);
        
        %%%%%%%%%%%%%%%%%
        % find unique vectors in the new labeling
        UQ1 = unique(labelVec);        
        for j = 1:numel(UQ1)
            % find the ith--jth subgroup
            idx1 = find((labelVec == UQ1(j)));
            % store label in new label matrix
            outLabel(idx0(idx1)) = curLabel;
            % increament the label count
            curLabel = curLabel + 1;
        end
        
    end
end

%{
%%%%%%%%%%%%%%%%%
% example with equals
func{1}.phi = @(x,para)label_1(x,para);
func{1}.para.value = 1;
data = [1 0 1 0 0 0 0 0 1]';
groups = ones(size(data));
groups = hLabel(data,groups,func);

func{1}.phi = @(x,para)label_2(x,para);
func{1}.para.value = 2;
func{2}.phi = @(x,para)label_2(x,para);
func{2}.para.value = 2;
func{3}.phi = @(x,para)label_2(x,para);
func{3}.para.value = 2;
func{4}.phi = @(x,para)label_2(x,para);
func{4}.para.value = 2;
groups = ones(size(STACK,2),1);
groups = hLabel(STACK',groups,func);



%}