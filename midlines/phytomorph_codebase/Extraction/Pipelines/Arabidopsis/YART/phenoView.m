function [] = phenoView(pheno,toLabel)    
    for e = 1:numel(pheno)
        vertLine = pheno(e).gamma(:,[1 end]);
        vertLine(2,2) = vertLine(2,1);
        
        horzLine = pheno(e).gamma(:,[1 end]);
        horzLine(1,2) = horzLine(1,1);
        
        plot(vertLine(2,:),vertLine(1,:),'m')
        plot(horzLine(2,:),horzLine(1,:),'m')
        
        plot(pheno(e).gamma(2,:),pheno(e).gamma(1,:),'r');
        plot(pheno(e).gamma(2,1),pheno(e).gamma(1,1),'go');
        plot(pheno(e).gamma(2,end),pheno(e).gamma(1,end),'go');
        plot(pheno(e).gamma(2,[1 end]),pheno(e).gamma(1,[1 end]),'b');
        quiver(pheno(e).gamma(2,end),pheno(e).gamma(1,end),pheno(e).E(2),pheno(e).E(1),20,'Color','c');
        %quiver(pheno(e).gamma(2,end),pheno(e).gamma(1,end),-pheno(e).E(2),-pheno(e).E(1),20,'Color','c');
        if toLabel
            loc = pheno(e).gamma(:,1);
            text(loc(2),loc(1),num2str(e),'BackgroundColor',[1 1 1]);
        end
        
        if ~isempty(pheno(e).LAT)
            for l = 1:numel(pheno(e).LAT)
                    phenoView(pheno(e).LAT(l),0);
            end
        end
    end
end