function [out1 out2] = testCondorFunction(arg1,arg2,arg3,arg4)
    out1 = arg1 + arg2;
    out2 = arg1 - arg2;
    fprintf(['Hello World.Recompile4.\n']);
    if nargin >= 3
        mkdir(arg3);
        image;
        saveas(gca,[arg3 'test.tif']);
    end
    if nargin == 4
        mkdir(arg4);
        image(zeros(10,10))
        saveas(gca,[arg4 'test2.tif']);
    end
    close all
end