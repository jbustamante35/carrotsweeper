% next attempt at prolog
%************************************
% make containers:contains:inverse
subClass(logicalContainer,container).
subClass(physicalContainer,container).
t(contains).
t(containedIn).
i(contains,containedIn).
contains(field,box).
contains(box,frogs).
