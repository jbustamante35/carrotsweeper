function [] = tensorView(data,para,h)
    try
        switch para.type
            case 'phytoApoint'     
                data = permute(data,[2 1 3]);
                for e = 1:size(data,3)
                    hp = plot(data(1,1,e),data(2,1,e),'.');
                    setProps(hp,para.props);
                end
            case 'phytoPoint'     
                data = permute(data,[2 1 3]);
                for e = 1:size(data,3)
                    hp = plot(data(1,3,e),data(2,3,e),'.');
                    setProps(hp,para.props);
                end
            case 'phytoAcurve'
                data = permute(data,[2 1 3]);
                for e = 1:size(data,3)
                    hp = plot(data(1,:,e),data(2,:,e));
                    setProps(hp,para.props);
                end
            case 'segmentSet'
                for e = 1:size(data,2)
                    vec = [data(1:2,e)';data(3:4,e)'];
                    hp = plot(vec(:,1),vec(:,2));
                    setProps(hp,para.props);
                end
            case 'phytoAaffine'
                for e = 1:size(data,3)
                    para.mag = 10;
                    pt = data(1:2,end,e);
                    for col = 1%:2%(size(data,1)-1)
                        d = data(1:2,col,e);
                        d = flipud(d);
                        if isfield(para,'mag')
                            d = d*para.mag;
                        end                    
                        vec = [pt';d'];
                        vec = cumsum(vec,1);                    
                        hp = plot(vec(:,1),vec(:,2),'b');
                        %setProps(hp,para.props);
                    end
                    hp = plot(pt(1),pt(2),'*');
                    setProps(hp,para.props);
                end
            case 'phytoAsequence'
                %%% curve
                curve = squeeze(data(:,:,end));
                curve = phytoAcurve(curve);
                curve.view(h,[]);
                for e = 1:size(data,1)
                    aSeq = squeeze(data(e,:,:));
                    if size(aSeq,2) >= 4
                        aSeq = [[aSeq(1:2,1:2) aSeq(1:2,end)];aSeq(end,[1 2 size(aSeq,2)])];
                    end
                    %R = [[cos(pi/2) sin(pi/2) 0];[-sin(pi/2) cos(pi/2) 0];[0 0 1]];
                    %aSeq = aSeq*R;
                    aSeq = phytoAaffine(aSeq);
                    aSeq.view(h,[]);
                end
            case 'phytoApatch'
                imshow(data);
            case 'phytoApatchSequence'
                for e = 1:size(data,1)
                    p = squeeze(data(e,:,:));
                    imshow(p,[]);
                    drawnow
                end
        end
    catch ME
        ME;
    end
end


function [] = setProps(h,props)
    try
        if ~isempty(props)
            set(h,props);
        end
    catch ME
        ME
    end
end