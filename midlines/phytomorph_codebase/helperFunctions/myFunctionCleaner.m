function [funcHandle] = myFunctionCleaner(funcHandle)
    s = functions(funcHandle);
    strF = func2str(funcHandle);
    
    
    
    f = fields(s.workspace{2});
    for blah = 1:numel(f)
        if isempty(strfind(strF,f{blah}))
            rm(blah) = 0;
        else
            rm(blah) = 1;
            CMD = [f{blah} '= s.workspace{2}.(f{blah});'];
            eval(CMD)
        end
    end
    f(logical(rm))
    clear f s ans strF blah CMD rm
    funcHandle =str2func(func2str(funcHandle));
    %{
    
    fidx1 = strfind(strF,'(');
    fidx2 = strfind(strF,')');
    args = strF((fidx1(end)+1):(fidx2(end)-1));
    %}
    
    
end