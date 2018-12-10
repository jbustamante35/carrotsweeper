function [] = mP(fibreBundle,figureHandle)
    defaultCL = 'b';
    defaultSZ =  6;
    %figure(figureHandle);
    for e = 1:numel(fibreBundle)
        L{e} = fibreBundle{e}.name;
        switch fibreBundle{e}.class
            case 'affine'
                cl = fibreBundle{e}.cl;
                sz = defaultSZ;
                alphaS = 50;
                plot(fibreBundle{e}.data(:,5),fibreBundle{e}.data(:,6),['k' 'o'],'MarkerSize',sz,'MarkerFaceColor',cl);
                quiver(fibreBundle{e}.data(:,5),fibreBundle{e}.data(:,6),alphaS*fibreBundle{e}.data(:,1),alphaS*fibreBundle{e}.data(:,2),0,'Color','g');
                quiver(fibreBundle{e}.data(:,5),fibreBundle{e}.data(:,6),alphaS*fibreBundle{e}.data(:,3),alphaS*fibreBundle{e}.data(:,4),0,'Color','b');
            case 'phytoApoint'
                cl = fibreBundle{e}.cl;
                sz = defaultSZ;
                plot(fibreBundle{e}.data(:,1),fibreBundle{e}.data(:,2),['k' 'o'],'MarkerSize',sz,'MarkerFaceColor',cl);
            case 'vectorField'
                cl = fibreBundle{e}.cl;
                sz = defaultSZ;
                alphaS = 1;
                plot(fibreBundle{e}.data(:,end-1),fibreBundle{e}.data(:,end),['k' 'o'],'MarkerSize',sz,'MarkerFaceColor',cl);
                quiver(fibreBundle{e}.data(:,end-1),fibreBundle{e}.data(:,end),alphaS*fibreBundle{e}.data(:,1),alphaS*fibreBundle{e}.data(:,2),0,'Color',cl);
            case 'phytoACurve'
                cl = fibreBundle{e}.cl;
                sz = defaultSZ;
                alphaS = 1;
                plot(fibreBundle{e}.data(:,1),fibreBundle{e}.data(:,2),'Color',cl);
        end
    end
    %legend(L);
end