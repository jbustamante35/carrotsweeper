function [germ frame] = score(datag,dataF,Ug,Eg,Uf,Ef,netG,netF,nf,hmm,thresholdV)


        fd = [];
        for e = 1:numel(dataF)
            WIN = im2col(dataF{e},[nf 1],'sliding');
            [C{e}] = PCA_REPROJ_T(WIN,Ef{e},Uf{e});
            fd = [fd;C{e}];
        end
        
        gd = [];
        for e = 1:numel(datag)
            [Cg{e}] = PCA_REPROJ_T(datag{e},Eg{e}(1:size(datag{e},1),:),Ug{e}(1:size(datag{e},1),:));
            gd = [gd;Cg{e}];
        end
        
        germ = sim(netG,gd);
        
        germ = germ(1,:) < thresholdV;
        %[~,germ] = max(germ,[],1);
        %germ = germ - 1;
        
        
        
        
        C{end+1} = sim(netF,fd);
        %C{end+1} = col2im(frameS,[nf 1],size(dataF{1}),'sliding');
        

        observation_labels = [];
        for e = 1:numel(C)
            observation_labels = [observation_labels;e*ones(size(C{e},1),1)];
        end


        
        
        

        SIG = {};
        for e = 1:numel(C)
            SIG{e} = [];
            for t = 1:size(C{e},1)
                tmp = col2im(C{e}(t,:),[nf 1],size(dataF{1}),'sliding');
                %tmp = imfilter(tmp,fspecial('average',[31 1]),'replicate');
                SIG{e} = cat(3,SIG{e},tmp);
            end
        end
        
        %frame = zeros(size(dataF{1}));
        parfor e = 1:size(SIG{1},2)
            tic
            STIM = [];
            for s = 1:numel(SIG)
                tmpV = squeeze(SIG{s}(:,e,:))';
                STIM = [STIM;tmpV];
            end
    
            gidx = hmm.Viterbi(STIM,observation_labels,1);
    
    
            
    
            fidx = find(gidx==2);
            
            tmp = zeros(size(dataF{1},1),1);
            if ~isempty(fidx)
                tmp(fidx:end) = 1;
            end
            frame{e} = tmp;
            toc/10*size(SIG{1},2)/60
        end
        frame = cell2mat(frame);
        %{
        earlyT = 20;
        for e = 1:size(frameS,2)
        
            sig = frameS(:,e);
            sig = sig / max(sig);
            sig = imfilter(sig,fspecial('average',[5 1]),'replicate');
            sig = sig / max(sig);
            osig = sig;
            pidx = find(sig == imdilate(sig,strel('disk',5)));
            sig(1:earlyT) = 0;
            ridx = sig(pidx) < .2;
            pidx(ridx) = [];
            frame(:,e) = zeros(size(frameS,1),1);
            if ~isempty(pidx)
                frame(pidx(1):end,e) = 1;
            end
        end
        %}
        %frame(:,germ==0) = 0;
end