function [error] = pp_ver0(imageFile)
    error = 0;
    try
        tm = clock;
        oPath = '/mnt/tetra/nate/chips2/';
        mkdir(oPath);
        % get the chip mask
        fprintf(['Start with mask \n']);
        [I B] = getChipMask_ver1(imageFile,5000,1000);
        B = logical(B);
        fprintf(['Done with mask \n']);

        %{
        [T D] = readtext('/home/nate/Downloads/R091027.txt');
        XYZ = [];
        for c = 15:302
            bv = ones(1,numel(T{c}));
            bv(strfind(T{c},' ')) = 0;
            R = regionprops(logical(bv),'PixelIdxList');
            tmp = [];
            for l = 1:3
                tmp = [tmp [str2num(T{c}(R(l+1).PixelIdxList))]];
            end
            XYZ = [XYZ;tmp];




        end

        colorP = I(403:1300,3321:5010,:);

        RGB = xyz2rgb(XYZ,'OutputType','uint8');
        SS = round([size(colorP,1) size(colorP,2)].*[1/12 1/22])*.95;
        SP = SS/2;
        [g1 g2] = ndgrid(linspace(SP(1),SP(1)+SS(1)*12,12),linspace(SP(2),SP(2)+SS(2)*22,22));
        %}



        % get the boundaries for the mask
        dB = bwboundaries(B,8,'noholes');
        R = regionprops(B,'PixelIdxList','BoundingBox');
        [FUN GRID IG IGM] = ringSample(I,B,dB,R,1500,0);
        %{
        imshow(I,[]);
        drawnow
        for e = 1:numel(R)
            rectangle('Position',R(e).BoundingBox,'EdgeColor','r');
        end
        %}


        [p,nm,ext] = fileparts(imageFile);
        for e = 1:numel(dB)
            saveFile = [oPath nm '--chip--' num2str(e) '.mat'];
            sam = FUN(:,:,:,e);
            grid = GRID(:,:,:,e);
            img = IG{e};
            msk = IGM{e};
            save(saveFile,'sam','grid','img','msk');
        end

        etm = etime(clock,tm);
        fprintf(['done image in:' num2str(etm) '\n']);
    catch ME
        ME
        error = 1;
    end
    
end