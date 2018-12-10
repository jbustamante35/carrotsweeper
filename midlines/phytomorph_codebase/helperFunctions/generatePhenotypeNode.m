function [graph] = generatePhenotypeNode(graph,data,name,linker)
    graph.(linker).name = name;
    graph.(linker).size = size(data);
    graph.(linker).data = data;
    graph.(linker).type = class(data);
end