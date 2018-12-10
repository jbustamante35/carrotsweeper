function [] = genQRsheets(csvFile)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the data for the QR codes
    d = readtext(csvFile);
    % get the rendering type
    toRender = d(1,:);
    % get the fields
    fields = d(2,:);
    % delete the first rows
    d(1:2,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create max number of rows and columns
    RMAX = 10;      % max rows
    CMAX = 3;       % max columns
    SKIPX = 2.63;   % skip in X
    SKIPY = 1;      % skip in Y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BORDER = [.5 .5];
    QRsize = [.9 .9];
    initPosition = BORDER + [0 (RMAX-1)*SKIPY];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init to max to start with reset
    CR = RMAX + 1;
    CC = CMAX + 1;
    
    
    
    
    % get the screen size
    SS = get(groot,'ScreenSize');
    % fraction the size
    SS = SS*.9;
    % set the width of the paper
    WID = SS(4)*8.5/11;
    % 
    PS = [WID SS(4)];
    FIGLOC = [0 0 PS];
   
   
    
    % init total QR count
    cnt = 1;
    toCo = {};
    TOTP = size(d,1);
    for e = 1:TOTP
        
        
        % if row count is greater than max rows - reset
        if CR > RMAX
            
            currentPos = initPosition;
            
            figSaved = 0;
            CR = 1;
            CC = 1;
            close all
            
            fig = figure;
            
            
           
            fig.PaperPosition = [0 0 8.5 11];
           
            fig.Position = FIGLOC;
            fig.Units = 'inches';
            fig.MenuBar = 'none';
            fig.DockControls = 'off';
            RATIO = [fig.Position(4)].*[11].^-1;
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create the string for the QR code
        string = '{';
        for f = 1:numel(fields)
            if ~isempty(strfind(toRender(f),'Q'))
                if ~ischar(d{e,f})
                    d{e,f} = num2str(d{e,f});
                end
                string = [string fields{f} '_' d{e,f} '}{'];
            end
        end
        string(end) = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create human string
        oneLine = true;
        Hstring = '';
        for f = 1:numel(fields)
            if ~isempty(strfind(toRender{f},'H'))
                if ~ischar(d{e,f})
                    d{e,f} = num2str(d{e,f});
                end
                if oneLine
                    Hstring = [Hstring {[fields{f} ':' strrep(d{e,f},'_','\_')]}];
                else
                    Hstring = [Hstring {[fields{f} ':']}  {strrep(d{e,f},'_','\_')}];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        tic
        qr = single(qrcode_gen(string,'ErrQuality','M','Version',1));
        %qr = single(qrcode_gen(string));
        %qr = single(qrcode_gen('a'));
        toc
        
        
        %imshow(qr,[])
        newAX = axes('Visible','off','Units','inches');
        img = image(cat(3,qr,qr,qr)*255,'Parent',newAX);
        set(newAX,'Visible','off');
        %set(img.Parent,'Xlim',[0 .75]);
        %set(img.Parent,'Ylim',[0 .75]);
        
        tx = text(size(qr,2)+2,.25,Hstring,'FontSize',8,'FontUnits','inch');
        tx.Units = 'inches';
        tx.Position = [[QRsize(1) + .25] .25]*RATIO;
        
        set(img.Parent,'Position',[currentPos.*RATIO QRsize.*RATIO]);
        %set(gca,'Visible','off');
        currentPos(1) = currentPos(1) + SKIPX;
        if CC == CMAX
            currentPos(2) = currentPos(2) - SKIPY;
            currentPos(1) = initPosition(1);
            CR = CR + 1;
            CC = 1;
        else
             CC = CC + 1;
        end
      
       
        string
        drawnow
        
        if mod(e,RMAX*CMAX) == 0
            toCo{end+1} = ['page' num2str(cnt) '.pdf'];
            print(gcf,'-dpdf',['page' num2str(cnt) '.pdf']);
            cnt = cnt + 1;
            close all
            figSaved = 1;
        end
    end
    % save last sheet if it is not full
    if ~figSaved
        toCo{end+1} = ['page' num2str(cnt) '.pdf'];
        print(gcf,'-dpdf',['page' num2str(cnt) '.pdf']);
        cnt = cnt + 1;
        close all
        figSaved = 1;
    end
    
    CMD = ['pdftk'];
    for e = 1:numel(toCo)
        CMD = [CMD ' '  toCo{e}];
    end
    CMD = [CMD ' cat output finalQR.pdf'];
    system(CMD),'-echo';
    drawnow
end

%{
   
    csvFile = '/home/nate/Documents/EDGAR.csv';
    csvFile = '/home/nate/Downloads/MBS 2016 QR code sheet-1.csv';
    csvFile = '/home/nate/Downloads/Mike.csv';
    csvFile = '/home/nate/Downloads/Anibus.csv';
    csvFile = '/home/nate/Downloads/density QR.full.csv';
    csvFile = '/home/nate/Documents/EDGAR.csv';
    %csvFile = '/home/nate/Downloads/qr_csv.csv - qr_csv.csv.csv';
    genQRsheets(csvFile);
%}