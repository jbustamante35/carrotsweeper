function genRoomSetup()


%{
    cmd = {'C:tare:S1' 'C:tare:S2' 'C:tare:S3' 'C:tare:S4'};
    for e = 1:numel(cmd)
        cmdI{e} = single(qrcode_gen(cmd{e},'ErrQuality','M','Version',1));
        imshow(imresize(cmdI{e},2),[])
        imwrite(double(imresize(cmdI{e},2,'nearest')),['pageCMD' num2str(e) '.tif']);
        %print(gcf,'-dpdf',['pageCMD' num2str(e) '.pdf']);
        close all
    end
    
%}
    
    

    step(1).qr = {'R:1' 'R:2' 'R:3' 'R:4'};
    
    step(2).qr = {'T:ear1' 'T:ear2' 'T:cob' 'T:kernel'};
    
    step(3).qr = {''};
    
    step(4).qr = {'S:1' 'S:2' 'S:3' 'S:4'};
    
    
    
    
    QRPos(1,:) = [4*72 1*72]
    QRPos(2,:) = [4*72 (1+2.5)*72]
    QRPos(4,:) = [4*72 (3.5+2.5+2.5)*72]

    
    
    cnt = 1;
    toCo = {};


    close all
    
    for e1 = 1:numel(step(1).qr)
        qr1 = single(qrcode_gen(step(1).qr{e1},'ErrQuality','M','Version',1));
        %qr1 = imresize(qr1,10,'nearest');
        for e2 = 1:numel(step(2).qr)
            qr2 = single(qrcode_gen(step(2).qr{e2},'ErrQuality','M','Version',1));
            %qr2 = interp2(qr2,2);
            for e3 = 1:numel(step(3).qr)
                for e4 = 1:numel(step(4).qr)
                    qr4 = single(qrcode_gen(step(4).qr{e4},'ErrQuality','M','Version',1));
                    %qr4 = interp2(qr4,2);
                    
                    
                    if e1 == e4
                        SS = get(groot,'ScreenSize');
                        SS = SS*.9;
                        WID = SS(4)*8.5/11;
                        PS = [WID SS(4)];
                        FIGLOC = [0 0 PS];

                        fig = figure;
                        figSaved = 0;
                        %fig.MenuBar = 'none';
                        %fig.DockControls = 'off';
                        %hold on
                        %u = fig.Units;
                        % fig.Units = 'inches';
                        fig.PaperPosition = [0 0 8.5 11];
                          fig.Position = FIGLOC;


                        toP = qr1;
                        newAX = axes('Visible','off','Units','pixels');
                        img = image(cat(3,toP,toP,toP)*255,'Parent',newAX);
                        set(newAX,'Visible','off');

                        set(img.Parent,'Xlim',[0 size(toP,2)]);
                        set(img.Parent,'Ylim',[0 size(toP,1)]);
                        set(img.Parent,'Position',[QRPos(1,:) [72 72]]);
                        text(size(qr1,2)-40,size(qr1,1)+5,'Step1','FontSize',20,'FontUnits','points','Rotation',90);
                        text(size(qr1,2)+20,size(qr1,1)+5,[step(1).qr{e1}] ,'FontSize',20,'FontUnits','points','Rotation',90);

                        toP = qr2;
                        newAX = axes('Visible','off','Units','pixels');
                        img = image(cat(3,toP,toP,toP)*255,'Parent',newAX);
                        set(newAX,'Visible','off');

                        set(img.Parent,'Xlim',[0 size(toP,2)]);
                        set(img.Parent,'Ylim',[0 size(toP,1)]);
                        set(img.Parent,'Position',[QRPos(2,:) [72 72]]);
                        text(size(qr1,2)-40,size(qr1,1)+5,'Step2','FontSize',20,'FontUnits','points','Rotation',90);
                        text(size(qr1,2)+20,size(qr1,1)+5,[step(2).qr{e2}] ,'FontSize',20,'FontUnits','points','Rotation',90);


                        toP = qr4;
                        newAX = axes('Visible','off','Units','pixels');
                        img = image(cat(3,toP,toP,toP)*255,'Parent',newAX);
                        set(newAX,'Visible','off');

                        set(img.Parent,'Xlim',[0 size(toP,2)]);
                        set(img.Parent,'Ylim',[0 size(toP,1)]);
                        set(img.Parent,'Position',[QRPos(4,:) [72 72]]);
                        text(size(qr1,2)-40,size(qr1,1)+5,'Step4','FontSize',20,'FontUnits','points','Rotation',90);
                        text(size(qr1,2)+20,size(qr1,1)+5,[step(4).qr{e4}] ,'FontSize',20,'FontUnits','points','Rotation',90);

                        drawnow
                        toCo{end+1} = ['page' num2str(cnt) '.pdf'];
                        print(gcf,'-dpdf',['page' num2str(cnt) '.pdf']);
                        cnt = cnt + 1;
                        close all
                    end
                end
            end
        end
    end
    
     
    CMD = ['pdftk'];
    for e = 1:numel(toCo)
        CMD = [CMD ' '  toCo{e}];
    end
    CMD = [CMD ' cat output finalQR_room.pdf'];
    system(CMD),'-echo';
    drawnow
        
    
end