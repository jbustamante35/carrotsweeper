function [f] = logicFunc(gate,data,target)
    tic
    gate = uint8(round(gate));
    res = any(bitand(data,data) > 0,2);
    fidx1 = find(target==1);
    fidx0 = find(target==0);
    rnd0 = randi(numel(fidx0),[15000 1]);
    f = -matthews_correlation(res([fidx1;fidx0(rnd0)]),target([fidx1;fidx0(rnd0)]));
    toc
end

%{

    close all
    h1 = figure;
    h2 = figure;
    TOT = [];
    POS = [];
    mC = [];
    gate = uint8(round(gate));
    cnt = 1;

    target1 = target;
    target1c = bitcmp(target1);
    %{
    m = 0;
    v = 2^m;
    sG = bitshift(target1,m);
    msk = repmat(uint8(v),[size(data,1) 1]);
    mskC = bitcmp(msk);
    Z = zeros(size(msk),'uint8');
    pos = 1;
    P = bitor(bitand(data(:,pos),Z),bitand(sG,msk));
    Q = bitor(bitand(data(:,pos),Z),bitand(data(:,pos),mskC));
    H = bitor(P,Q);
    data(:,pos) = H;
    %}
    %basis = zeros(size(data,1),1,'logical');
    basis = [];
    disp = true;
    c = [];
    gc = [];
    for pos = 1:size(data,2)
        for m = 0:7

            tic;

            v = 2^m;
            
            %if pos == 1
                %tmpTB(:,m+1) = bitshift(target1,m);
            %end
            %tmpTc = bitshift(target1c,m);

            tmpT = bitshift(target1,m);
            %tmpT = tmpTB(:,m+1);

            msk = repmat(uint8(v),[size(data,1) 1]);

            
            F = bitand(msk,data(:,pos));

            
            F1_bin = F == v;
            F0_bin = ~F1_bin;


            F1_bin = uint8(F1_bin);
            F0_bin = uint8(F0_bin);



            tmpT1_bin = uint8(tmpT == v);
            tmpT1c_bin = uint8(~tmpT1_bin);


            %[EN(cnt)] = conditionalEntropy(tmpT1_bin,F1_bin);
            %[ME(cnt)] = mutualEntropy(F1_bin,tmpT1_bin);
            ME(cnt) = mutualinfo(F1_bin,tmpT1_bin);

            %[c(cnt),C] = informationContent(tmpT1_bin,F1_bin);
            

            if ~isempty(basis)
                %[gc(cnt),gC,totC(cnt)] = informationContentGain(tmpT1_bin,F1_bin,basis);

            end




            %MCCS(:,:,cnt) = C;




            cnt = cnt + 1;
            if disp
                figure(h1)
                plot(ME)


%{
                figure(h2)
                plot(c,'g');
                hold on
                plot(gc,'b');
                figure(h1)
                hold on
                %plot(totC,'r');
                hold off
%}
                drawnow
            end
            


            cnt


            toc
        end
    end
    




    [scC,sIDX] = sort(ME,'descend');
    [b,c,r] = ndgrid(1:8,1:13,1:13);
    px = [r(:) c(:) b(:)];

    S = zeros(size(target),'uint8');
    for s = 1
        m = px(sIDX(s),3) - 1;
        pos = (px(sIDX(s),2)-1)*8 + 4;
        v = 2^m;
        msk = repmat(uint8(v),[size(data,1) 1]);
        F = bitand(msk,data(:,pos));
        F = bitshift(F,-(m)+(s-1));
        S = bitor(S,F);
    end
    basis = S==1;


    [scC,sIDX] = sort(c,'descend');
    [b,c,r] = ndgrid(1:8,1:13,1:13);
    px = [r(:) c(:) b(:)];


    S = zeros(size(target),'uint8');
    for s = 1
        m = px(sIDX(s),3) - 1;
        pos = (px(sIDX(s),2)-1)*8 + 4;
        v = 2^m;
        msk = repmat(uint8(v),[size(data,1) 1]);
        F = bitand(msk,data(:,pos));
        F = bitshift(F,-(m)+(s-1));
        S = bitor(S,F);
    end
    basis = S==1;




    [scC,sIDX] = sort(gc,'descend');
    [b,c,r] = ndgrid(1:8,1:13,1:13);
    px = [r(:) c(:) b(:)];

    for s = 1
        m = px(sIDX(s),3) - 1;
        pos = (px(sIDX(s),2)-1)*8 + 4;
        v = 2^m;
        msk = repmat(uint8(v),[size(data,1) 1]);
        F = bitand(msk,data(:,pos));
        F = bitshift(F,-(m)+(s-1));
        S = bitor(S,F);
    end
    basis = [basis S==1];



    
    pos = floor(double(S)*8^-1)+1;
    m = rem(double(S)-1,8)+1;
    v = 2.^m;
    muxBank = zeros(size(data,1),32,'uint8');
    mIDX = sub2ind(size(muxBank),(1:size(muxBank,1))',pos);
    muxBank(mIDX) = v;


            
    F1_bin = S == 255;
    F0_bin = ~F1_bin;


    F1_bin = single(F1_bin);
    F0_bin = single(F0_bin);

    tmpT1_bin = single(target == 1);
    tmpT1c_bin = single(~tmpT1_bin);



    TP = mtimesx(tmpT1_bin,'T',F1_bin);
    FP = mtimesx(tmpT1c_bin,'T',F1_bin);

    TN = mtimesx(tmpT1c_bin,'T',F0_bin);
    FN = mtimesx(tmpT1_bin,'T',F0_bin);


    C(1,1) = TN;
    C(1,2) = FN;
    C(2,1) = FP;
    C(2,2) = TP;




    TOTP = sum(target);
    TOTN = sum(target==0);
    TOTPN = size(target,1);

tmp = C(1,:);



    for r = 1:size(C,1)
        tmp = bitand(C(r,:),tmp);
        r
    end
%}
