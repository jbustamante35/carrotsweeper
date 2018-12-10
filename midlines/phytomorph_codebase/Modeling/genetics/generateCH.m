function [C] = generateCH(cL,maxS,type)
    C.ch = [];
    C.cn = [];
    for e = 1:numel(cL)
        switch type
            case 'binary'
                 tmp = [zeros(cL(e),1);ones(cL(e),1)];
            case 'inbred'
                tmp = randi(maxS,cL(e),1);
                tmp = [tmp;tmp];
            case 'rand'
                tmp = [randi(maxS,cL(e),1);randi(maxS,cL(e),1)];
        end
        ibd = ((.5-eps)).*rand(size(tmp));
        
        C.ch = [C.ch;tmp+ibd];
        C.cn = [C.cn;[e*ones(2*cL(e),1) [zeros(cL(e),1);ones(cL(e),1)]]];
    end
end