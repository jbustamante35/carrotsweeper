function [diameter] = measureStemDiameter(MASK,STEM_SNIP,outLN,basePoint)
    if ~isempty(basePoint)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the stem diameter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        STEM = MASK((end-STEM_SNIP):end,:);
        STEM = logical(STEM) - ~imfill(~logical(STEM),[size(STEM,1) round(basePoint(2))]);
        diameter = sum(STEM,2);
        [diameter,sidx] = sort(diameter);
        %{
        sdiameter = std(diameter(outLN:end-(outLN-1)),1,1);
        diameter = mean(diameter(outLN:end-(outLN-1)));
        %}
        PER = 2;
        MD = mode(diameter);
        diameter(diameter > MD*PER) = [];
        diameter = mean(diameter);
        %{
        if sdiameter < 1;sdiameter=1;end
        
        fdiameter = mean(diameter((end-outLN):end));
        diameter = mean(diameter(outLN:end-(outLN-1)));
        
        
        if (fdiameter - diameter)/sdiameter > 3 || sdiameter > 5
            diameter = inf;
        end
        %}
        
        %{
        if any(sidx((end-outLN):end) > (size(STEM,1) - outLN/2))
            diameter = inf;
        else
            diameter = mean(diameter(outLN:end-(outLN-1)));
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the stem diamter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        diameter = NaN;
    end
end