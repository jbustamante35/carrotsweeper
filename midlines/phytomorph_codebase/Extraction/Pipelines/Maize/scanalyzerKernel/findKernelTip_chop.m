function [mm1,mm2,mm3,idx1,idx2] = findKernelTip_chop(dB,L,N,I,disp)













    try

        
        
        if disp
            h1 = figure;
            h2 = figure;
        end

        a = zeros(size(dB,1),1);

        for e = 1:size(dB,1)
            
            
            [TAN,NOR,MID,seg,mm1,mm2,mm3,idx1,idx2] = getKframe(e,L,dB);
            [TAN1,NOR1,MID1,seg1,mm11,mm21,mm31,idx11,idx21] = getKframe(e-N,L,dB);
            [TAN2,NOR2,MID2,seg2,mm12,mm22,mm32,idx12,idx22] = getKframe(e+N,L,dB);
            
            delta1 = MID1 - MID;
            delta2 = MID2 - MID;
            
            MM31 = bsxfun(@minus,mm31,delta1);
            MM32 = bsxfun(@minus,mm32,delta2);
            MM21 = bsxfun(@minus,mm21,delta1);
            MM22 = bsxfun(@minus,mm22,delta2);
            
            MM22 = diff(MM22,1,1);
            MM21 = diff(MM21,1,1);
            MM2 = diff(mm2,1,1);
            MM22 = MM22 / norm(MM22);
            MM21 = MM21 / norm(MM21);
            MM2 = MM2 / norm(MM2);
            a1 = acos(MM2*MM21')*180/pi;
            a2 = acos(MM2*MM22')*180/pi;
            a(e) = mean([a1 a2]);
            IDX(e) = idx2;
            %{
            if disp
                figure(h1);
                subplot(2,1,1);
                imshow(I,[]);
                hold on
                plot(dB(:,2),dB(:,1),'m');
                hold on
                plot(seg(:,2),seg(:,1),'r');
                %plot(seg2(:,2),seg2(:,1),'w');
                plot(mm1(:,2),mm1(:,1),'k');
                plot(mm2(:,2),mm2(:,1),'b');
                plot(mm3(:,2),mm3(:,1),'c');
                plot(mm2(:,2),mm2(:,1),'b');
                plot(mm3(:,2),mm3(:,1),'c');
                plot(dB(idx1,2),dB(idx1,1),'b*');
                plot(dB(idx2,2),dB(idx2,1),'c*');


                plot(mm31(:,2),mm31(:,1),'c');
                plot(mm32(:,2),mm32(:,1),'c');

                plot(mm21(:,2),mm21(:,1),'b');
                plot(mm22(:,2),mm22(:,1),'b');
                hold off
                axis equal
                
                figure(h1);
                subplot(2,1,2)
                plot(a);
                drawnow
            end
            %}
        end
        
        [~,midx] = max(a);
        
        
            
        [TAN,NOR,MID,seg,mm1,mm2,mm3,idx1,idx2] = getKframe(midx,L,dB);
    catch ME
        ME
    end

    
    
    
    
end



function [TAN,NOR,MID,seg,mm1,mm2,mm3,idx1,idx2] = getKframe(e,L,dB)

   

    str = e;
    
    if str < 1
        str = size(dB,1) + str;
    end
    
    stp = str + L;
    flag = false;
    if stp >= size(dB,1)
        flag = true;
    end
    stp = mod(stp,size(dB,1))+1;
    if ~flag
        seg = dB(str:stp,:);
    else
        seg = [dB(str:end,:);dB(1:stp,:)];
    end
    TAN = seg(end,:) - seg(1,:);
    NOR = [TAN(2) -TAN(1)];
    MID = seg(1,:) + .5*TAN;
    NOR = NOR / norm(NOR);


    mag = 1600;
    mm1 = [seg(1,:);seg(end,:)];
    mm2 = [MID;MID+mag*NOR];
    mm3 = [MID;MID-mag*NOR];

    l2 = [linspace(mm2(1,1),mm2(2,1),1000)' linspace(mm2(1,2),mm2(2,2),1000)'];
    l3 = [linspace(mm3(1,1),mm3(2,1),1000)' linspace(mm3(1,2),mm3(2,2),1000)'];
    [D2,I2] = pdist2(dB,l2,'euclidean','Smallest',1);
    [D3,I3] = pdist2(dB,l3,'euclidean','Smallest',1);
    [~,s2] = min(D2);
    [~,s3] = min(D3);

    mm2 = [MID;dB(I2(s2),:)];
    mm3 = [MID;dB(I3(s3),:)];

    idx1 = I2(s2);
    idx2 = I3(s3);
    
    
end