function [G] = issueGrades(data,iIDX,iE,iU,iPDF_U,iPDF_C,eIDX,eE,eU,ePDF_U,ePDF_C,rEI,rIE,idU,idE,myiPDF,myePDF)
    for e = 1:size(data,1)
        tic
        tmpdata = data(e,:);
        tmpdata = zscore(tmpdata,1,2);
        %tmpdata = reshape(tmpdata,[50 50]);
        %tmpdata = mean(cat(3,tmpdata,flipdim(tmpdata,1),flipdim(tmpdata,2),flipdim(flipdim(tmpdata,1),2)),3);
        
       
            tic
            itmp = tmpdata(iIDX);
            itmp = itmp - mean(itmp);
            iC = PCA_REPROJ(itmp,iE,iU);
            iS = PCA_BKPROJ(iC,iE,iU);

            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            iDelta = (iS - itmp);
            iCD = PCA_REPROJ(iDelta,idE,idU);
            iGD = [iCD];
            if ~isempty(myiPDF)
                tmpd = [];
                for d = 1:numel(iGD)
                    tmpd(d) = interp1(myiPDF.XD(d,:),myiPDF.FD(d,:),iGD(d),'linear',0);
                end
                iGPD(e) = prod(tmpd);
            else
                iGPD(e) = mvnpdf(iG,iPDF_U,iPDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}

            iEr = norm(itmp - iS);
            iG = [iC iEr];
            iGT = myiPDF.MAP1*[iC(1:2) ones(size(iG,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myiPDF)
                tmp = [];
                for d = 1:numel(iG)
                    tmp(d) = interp1(myiPDF.X(d,:),myiPDF.F(d,:),iG(d),'linear',0);
                end
                tmp2 = ba_interp2(myiPDF.X2Yi,iGT(2),iGT(1),'linear');
                iGP2(e) = prod(tmp2)*tmp(end);
                iGP(e) = prod(tmp);
            else
                iGP(e) = mvnpdf(iG,iPDF_U,iPDF_C);
            end
            toc
            %{
            
            f(e,:) = iGT;
            a(e,:) = tmp;
            b(e,:) = iG;
            c(e,:) = tmp2;
            g(e,:) = tmpd;
            h(e) = iU*itmp'/norm(iU);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            freeWheel = zeros(size(tmpdata));
            freeWheel(iIDX) = iS;

            freeWheelD = zeros(size(tmpdata));
            freeWheelD(iIDX) = iDelta;
            %freeWheel(eIDX) = eS;
            W(:,:,e) = reshape(freeWheel,[50 50]);
            WD(:,:,e) = reshape(freeWheelD,[50 50]);
            sI(:,:,e) = reshape(tmpdata,[50 50]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %}
            tic
            etmp = tmpdata(eIDX);
            eC = PCA_REPROJ(etmp,eE,eU);
            eS = PCA_BKPROJ(eC,eE,eU);
            eEr = sum((eS - etmp).^2,2).^.5;
            eG = [eC eEr];
            eGT = myePDF.MAP2*[eC(1:2) ones(size(eC,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myePDF)
                tmp = [];
                for d = 1:numel(eG)
                    tmp(d) = interp1(myePDF.X(d,:),myePDF.F(d,:),eG(d),'linear',0);
                end
                tmp2 = ba_interp2(myePDF.X2Ye,eGT(2),eGT(1),'linear');
                eGP2(e) = prod(tmp2)*tmp(end);
                eGP(e) = prod(tmp);
            else
                eGP(e) = mvnpdf(eG,ePDF_U,ePDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc



            tic
            epC = (iC*rIE);
            epS = PCA_BKPROJ(epC,eE,eU);
            epEr = sum((epS - etmp).^2,2).^.5;
            epG = [epC epEr];
            epGT = myePDF.MAP2*[epC(1:2) ones(size(epC,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myePDF)
                tmp = [];
                for d = 1:numel(epG)
                    tmp(d) = interp1(myePDF.X(d,:),myePDF.F(d,:),epG(d),'linear',0);
                end
                tmp2 = ba_interp2(myePDF.X2Ye,epGT(2),epGT(1),'linear');
                epGP2(e) = prod(tmp2)*tmp(end);
                epGP(e) = prod(tmp);
            else
                epGP(e) = mvnpdf(eG,ePDF_U,ePDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc


            tic
            ipC = (eC*rEI);
            ipS = PCA_BKPROJ(ipC,iE,iU);
            ipEr = sum((ipS - itmp).^2,2).^.5;
            ipG = [ipC ipEr];
            ipGT = myiPDF.MAP1*[ipC(1:2) ones(size(ipC,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myiPDF)
                tmp = [];
                for d = 1:numel(ipG)
                    tmp(d) = interp1(myiPDF.X(d,:),myiPDF.F(d,:),ipG(d),'linear',0);
                end
                tmp2 = ba_interp2(myiPDF.X2Yi,ipGT(2),ipGT(1),'linear');
                ipGP2(e) = prod(tmp2)*tmp(end);
                ipGP(e) = prod(tmp);
            else
                ipGP(e) = mvnpdf(ipG,iPDF_U,iPDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc


            tic
            for itr = 1:15
                eC = .5*(epC + eC);
                ipC = (eC*rEI);
                iC = .5*(ipC + iC);
                epC = (iC*rIE);
                %{
                store(itr,:) = iC;
                ifS = PCA_BKPROJ(iC,iE,iU);
                efS = PCA_BKPROJ(eC,eE,eU);
                ifEr = sum((ifS - itmp).^2,2).^.5;
                efEr = sum((efS - etmp).^2,2).^.5;
                ifG = [iC ifEr];
                ifGP(itr) = mvnpdf(ifG,iPDF_U,iPDF_C);
                efG = [eC efEr];
                efGP(itr) = mvnpdf(efG,ePDF_U,ePDF_C);
                %}
            end
            toc



            tic
            ifS = PCA_BKPROJ(iC,iE,iU);
            efS = PCA_BKPROJ(eC,eE,eU);
            ifEr = sum((ifS - itmp).^2,2).^.5;
            efEr = sum((efS - etmp).^2,2).^.5;
            ifG = [iC ifEr];
            ifGT = myiPDF.MAP1*[iC(1:2) ones(size(iC,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myiPDF)
                tmp = [];
                for d = 1:numel(ifG)
                    tmp(d) = interp1(myiPDF.X(d,:),myiPDF.F(d,:),ifG(d),'linear',0);
                end
                tmp2 = ba_interp2(myiPDF.X2Yi,ifGT(2),ifGT(1),'linear');
                ifGP2(e) = prod(tmp2)*tmp(end);
                ifGP(e) = prod(tmp);
            else
                ifGP(e) = mvnpdf(ifG,iPDF_U,iPDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc



            tic
            efG = [eC efEr];
            efGT = myePDF.MAP2*[eC(1:2) ones(size(eC,1),1)]';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~isempty(myePDF)
                tmp = [];
                for d = 1:numel(efG)
                    tmp(d) = interp1(myePDF.X(d,:),myePDF.F(d,:),efG(d),'linear',0);
                end
                tmp2 = ba_interp2(myePDF.X2Ye,efGT(2),efGT(1),'linear');
                efGP2(e) = prod(tmp2)*tmp(end);
                efGP(e) = prod(tmp);
            else
                efGP(e) = mvnpdf(efG,ePDF_U,ePDF_C);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            toc

            %{
            fM = zeros(50,50);
            fM(iIDX) = ifS;
            fM(eIDX) = efS;
            imshow(fM);
            %}

        
        toc*(size(data,1)-e)/60/10
    end
    G = [iGP' eGP' ipGP' epGP' ifGP' efGP' iGP2' eGP2' ipGP2' epGP2' ifGP2' efGP2'];
   
end