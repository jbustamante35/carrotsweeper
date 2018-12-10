function [hmm] = makeChainRepeat_foEmergence(U0,C0,U1,C1,D,HARDLINE_HOLD,n,nCOMP,gmmNUM)


    hmm = my_hmm();
    hmm.dn = 1*ones(nCOMP,1);
    
    
    for u = 1:n
        
        
        Dis0{u} = myProb(U0',C0);
        Dis1{u} = myProb(U1',C1);
        
        %{
        Dis0{u}.fitToKS(D{1}',[0 1],256);
        Dis1{u}.fitToKS(D{2}',[0 1],256);
        %}
        %{
        if gmmNUM(1) > 0
            Dis0{u}.fitToGMM(D{1}',2,1);
        end
        
        if gmmNUM(2) > 0
            Dis1{u}.fitToGMM(D{2}',2,1);
        end
        %}
        
        
        n0{u} = hmm_node(['backGround-' num2str(u)]);
        n1{u} = hmm_node(['foreGround-' num2str(u)]);
        
        var00{u} = constantTransitionFunction(HARDLINE_HOLD);
        var01{u} = constantTransitionFunction(1-HARDLINE_HOLD);
        
        var10{u} = constantTransitionFunction(1-HARDLINE_HOLD);
        
        %var11{u} = VarheavisideTransitionFunction(800,@(x,y)lt(x,y),HARDLINE_HOLD);
        var11{u} = constantTransitionFunction(HARDLINE_HOLD);
        
        n0{u}.attachDistribution(Dis0{u},1);
        n1{u}.attachDistribution(Dis1{u},1);
        
        n0{u}.attachNode(n1{u},var01{u});
        n0{u}.attachNode(n0{u},var00{u});
        
        n1{u}.attachNode(n1{u},var11{u});
        
        if u ~= 1
            n1{u-1}.attachNode(n0{u},var10{u});
        end
        
       
        hmm.addNode(n0{u});
        hmm.addNode(n1{u});
       
    end
    u = u + 1;
    n0{u} = hmm_node(['backGround-' num2str(u)]);
    Dis0{u} = myProb(U0',C0);
    var00{u} = constantTransitionFunction(1);
    var10{u} = constantTransitionFunction(1-HARDLINE_HOLD);
    n0{u}.attachDistribution(Dis0{u},1);
    n0{u}.attachNode(n0{u},var00{u}); 
    n1{u-1}.attachNode(n0{u},var10{u});
    hmm.addNode(n0{u});
    
end