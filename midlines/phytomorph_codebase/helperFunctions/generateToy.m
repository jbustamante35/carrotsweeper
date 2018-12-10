function [X,Y,E] = generateToy(type,para)

    switch type
        case 'linear'
            X = pullRandom('uniform',{20 40},sz)
            
            
            N = para.N;
            mb = rand(para.ndims,2);
            X = rand(N,para.ndims);
    end

end

function [r] = pullRandom(distribution,para,sz)
    r = random(distribution,para{:});
end

%{
    [X,Y,E] = generateToy('linear',para);
%}