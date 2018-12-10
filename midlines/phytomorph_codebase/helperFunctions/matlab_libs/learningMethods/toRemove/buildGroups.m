function [grpStruct] = buildGroups(data,para)
    
    %%%%%%%%%%%%%%%%%%%%
    % here is a list of methods for grouping the data
    %%%%%%%%%%
    % group on codomain - use recursive hLabel
    para.value = 1;
    func{1}.phi = @(x)label_1(x,para);
    % call recursive grouping
    groups = ones(size(data.codomain));
    groups = hLabel(data.codomain,groups,func);    
    grpStruct(1).groups = groups;
    grpStruct(1).name = 'group on codomain';
    
    
    
    %%%%%%%%%%%%%%%%%%%%
    % group on codomain - use recursive hLabel
    para.value = 3;
    func{1}.phi = @(x)kmeans(x,para.value);
    func{2}.phi = @(x)kmeans(x,para.value);
    % call recursive grouping
    groups = ones(size(data.domain,1),1);
    groups = hLabel(data.domain,groups,func);    
    
    grpStruct(1).groups = groups;
    grpStruct(1).name = 'group on codomain';
end