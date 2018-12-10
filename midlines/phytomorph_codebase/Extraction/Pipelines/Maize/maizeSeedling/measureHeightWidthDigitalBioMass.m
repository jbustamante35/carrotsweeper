function [plantHEIGHT,HEIGHT, WIDTH,dBIOMASS,hMassDistribution,vMassDistribution,CenterOfMass,StdDis] = measureHeightWidthDigitalBioMass(MASK)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure and spool height, width digitalBioMass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['starting: measure height,width,dBioMass,etc...\n']);
    
    %%% make Height measurement
    % sum the mask for the height calculation
    hMassDistribution = sum(MASK,2);
    % find the pixels
    fidx = find(hMassDistribution);
    if isempty(fidx)
        fidx = size(MASK,1);
    end
    % find the top pixel
    HEIGHT = fidx(1);
    
    %%% make Height measurement
    
    
    
    %%% make Width measurement
    vMassDistribution = sum(MASK,1);
    fidx = find(vMassDistribution);
    if isempty(fidx)
        fidx = 0;
    end
    WIDTH = max(fidx) - min(fidx);
    %%% make Width measurement
    
    
    
    
    %%% others
    % plant height
    plantHEIGHT = size(MASK,1) - HEIGHT;
    % find the biomass
    dBIOMASS = sum(MASK(:));
    [x,y] = find(MASK);
    [S C CenterOfMass E L ERR LAM] = PCA_FIT_FULL([x,y],2);
    [C] = PCA_REPROJ([x y],eye(2),CenterOfMass);
    StdDis = std(C,1,1);
    %%% others
    
    
    fprintf(['ending: measure height,width,dBioMass,etc...\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure and spool height, width digitalBioMass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end