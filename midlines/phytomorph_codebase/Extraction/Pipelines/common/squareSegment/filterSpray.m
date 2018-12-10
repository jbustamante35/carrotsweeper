function [filteredSpray] = filterSpray(X,filterValue)
    filteredSpray = X{1}(:,X{2}==filterValue);
end