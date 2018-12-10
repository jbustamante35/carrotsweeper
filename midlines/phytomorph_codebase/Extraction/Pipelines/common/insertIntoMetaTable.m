function [table,frameKey,uniqueKey] = insertIntoMetaTable(table,data,type,frameKey)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if key word image is present in the type
    % then create a frame for the image
    if contains(lower(type),'image')
        [table,frameKey] = insertIntoMetaTable(table,eye(3),'iFrame',frameKey);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if no frameKey is given, then create a random frame key
    % this should only be called by the above when an image is inserted
    % therefore a frame is needed
    if isempty(frameKey)
        %%%%%%%%%%%%%%%%%%%%%%%%
        % generate random key
        b = 2;
        R = randi(10^b,1);
        T = datenum(now);
        p = log10(T);
        p = floor(p);
        p = p + b;
        kR = (10^(p-1))*R;
        frameKey = kR + T;
        fprintf(['fk:' num2str(frameKey) '=' num2str(T) '+' num2str(kR) '\n']);
        %%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate random index key
    b = 4;
    R = randi((2^32)-1,1,'uint32');
    T = datenum(datetime);
    %{
    T = datenum(datetime);
    p = log10(T);
    p = floor(p);
    p = p + b;
    kR = (10^(p-1))*R;
    uniqueKey = kR + T;
    %}
    uniqueKey = R;
    %fprintf(['ik:' num2str(uniqueKey) '=' num2str(T) '+' num2str(kR) '\n']);
    fprintf(['ik:' num2str(uniqueKey) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if isempty(frameKey)
        frameKey = uniqueKey;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % insert into the table
    pr = size(table,1) + 1;
    table.key(pr) = uniqueKey;
    table.time(pr) = T;
    table.data{pr} = data;
    table.type{pr} = type;
    table.frameKey(pr) = frameKey;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end