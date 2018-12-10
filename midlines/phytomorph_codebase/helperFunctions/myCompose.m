function [fout] = myCompose(func1,func2,arg1,arg2,newArg2)


    in_1 = argnames(func1);
    outArgs = setdiff(in_1,arg1,'Stable');
    
    
    in_2 = argnames(func2);
    in_2 = subs(in_2,arg2,newArg2);
    
    
    in_2 = setdiff(in_2,outArgs,'Stable');
    
    outArgs = [outArgs(:);in_2(:)];
    
    fout(outArgs(:)) = compose(func1,func2,arg1,arg2,newArg2);
    %fout(outArgs(:)) = subs(func1,arg1,func2);
end