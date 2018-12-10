function [cmd] = gen_myf(myf)    
    for f = 1:numel(myf)
        inSTR = '(';
        for i = 1:size(myf(f).arg,2)
            inSTR = [inSTR 'myf(' num2str(f) ').arg{' num2str(i) '},'];        
        end
        inSTR(end) = [];
        inSTR = [inSTR ')'];    
        cmd{f} = [myf(f).f_handle inSTR ';'];
    end
end