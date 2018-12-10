function [U,V,D] = diskPCA(fileList,toVecFunc,selVecFunc,n)
    % get the mean and cov for processing
    [U,C] = diskMEANandCOV(fileList,toVecFunc,selVecFunc);
    [V,D] = eigs(C,n);
end

%{
    % generate sim set
    func = @(X)func_depthStack(X,0);
    tester = [];
    for e = 1:10
        fileList{e} = rand(200,200,3);
        tester = [tester func(fileList{e})];
    end
    [U,V,D] = diskPCA(fileList,func,3);
    [~,~,Ut,Et,] = PCA_FIT_FULL_T(tester,3);

    
%}