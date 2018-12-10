function [] = genQRsheets_latex(csvFile,oPath)

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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % init to max to start with reset
    CR = RMAX + 1;
    CC = CMAX + 1;
    
    
    headerString = ['\\documentclass{article}\n'...
                    '\\usepackage{amsmath}\n'...
                    '\\usepackage{graphicx}\n'...
                    '\\usepackage{xparse}\n'...
                    '\\usepackage{microtype}\n'...
                    '\\usepackage[margin=.75cm,showframe=false]{geometry}\n'...
                    '\n'...
                    '\\begin{document}\n'];
                    %'\\FloatBarrier\n'];
    
    pageHeaderString = ['\\begin{table}[t]\n'...
                    '\\begin{center}\n'...
                    '\\begin{tabular}{ ccc }\n'];
   
    pageFooterString = ['\\end{tabular}\n'...
                    '\\end{center}\n'...
                    '\\end{table}\n'];
    
                    
    %{  
    pageHeaderString = ['\\begin{center}\n'...
                    '\\begin{tabular}{ ccc }\n'];
   
    pageFooterString = ['\\end{tabular}\n'...
                   '\\end{center}\n'];
    %}
                    
               
    documentFooterString = ['\\end{document}\n'];
                
    qrString = ['\\begin{minipage}[c][1in][c]{1in}\n'...
               '\\centering \\raisebox{-.65in\\relax}{\\includegraphics[width=.9in,height=.9in]{#IMAGENAME#}}\n'...
               '\\end{minipage}\n'...
               '\\begin{minipage}[c][1in][c]{1.625in}\n'...
               '#HUMANTEXT# \n'...
                '\\end{minipage}\n'];
    
    % init total QR count
    cnt = 1;
    toCo = {};
    TOT = size(d,1);
    
    mkdir(oPath);
    
    pageNumber = 0;
    labelNumber = 0;
    
    pageLATEXfile = [oPath 'masterDoc.tex'];
    fid = fopen(pageLATEXfile,'wt');
    fprintf(fid,headerString);
    
    
    PART = RMAX*CMAX - rem(TOT,RMAX*CMAX);
    
    %PART = PART*double((floor(TOT*(RMAX*CMAX)^-1) >= 1));
    PART = PART*double((RMAX*CMAX ~= PART));
    
    FULLTOT = TOT + PART;
    
    fprintf(['Generating :Total:Found:Part:@:' num2str(FULLTOT) ':' num2str(TOT) ':' num2str(PART) ':\n'])
    
    for e = 1:FULLTOT
        fprintf(['start generating QR label latex code:' num2str(e) '\n']);
        
        % if row count is greater than max rows - reset
        if CR > RMAX
            figSaved = 0;
            fprintf(fid,pageHeaderString);
            pageNumber = pageNumber + 1;
            pagePath = [oPath];
            mkdir(pagePath);
            CR = 1;
            CC = 1;
        end

        if e <= TOT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % create the string for the QR code
            string = '{';
            for f = 1:numel(fields)
                if contains(toRender(f),'Q')
                    if ~ischar(d{e,f})
                        d{e,f} = num2str(d{e,f});
                    end
                    string = [string fields{f} '_' d{e,f} '}{'];
                end
            end
            string(end) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            string = '{noKEY_noVALUE}';
        end

        if e <= TOT
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
                        Hstring = [Hstring [fields{f} ': \newline ' strrep(d{e,f},'_','\_')] ' \newline '];
                    else
                        Hstring = [Hstring {[fields{f} ':']}  {strrep(d{e,f},'_','\_')}];
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            Hstring = 'Blank Label';
        end
        
        
        
        tic
        qr = single(qrcode_gen(string,'ErrQuality','M','Version',1));
        qr = interp2(qr,4,'nearest');
        tmpQRfile = [pagePath num2str(labelNumber) '.png'];
        imwrite(qr,tmpQRfile);
        
        
        
        TMPqrString = strrep(qrString,'#IMAGENAME#',[pagePath num2str(labelNumber) '.png']);
        
        Hstring = strrep(Hstring,'#','\#');
        Hstring = strrep(Hstring,'&','\&');
        Hstring = strrep(Hstring,'%','\%');
        Hstring = strrep(Hstring,'$','\$');
        Hstring = strrep(Hstring,'_','\_');
        Hstring = strrep(Hstring,'{','\{');
        Hstring = strrep(Hstring,'}','\}');
        Hstring = strrep(Hstring,'~','\~');
        Hstring = strrep(Hstring,'^','\^');
        Hstring = strrep(Hstring,'\','\\');
        
        
        %Hstring = 'helloWorld';
        TMPqrString = strrep(TMPqrString,'#HUMANTEXT#',Hstring);
        fprintf(fid,TMPqrString);
        
        
       
        labelNumber = labelNumber + 1;
        
        
        % handle the number of columns
        if CC == CMAX
            % drop down a row with the \\\ string
            fprintf(fid,'\\\\\n');
            % increment the row count
            CR = CR + 1;
            % reset the column count
            CC = 1;
            % 
        else
            % increment the column count
            %fprintf(fid,'\\end{minipage} \n & \n');
            fprintf(fid,'&\n');
            CC = CC + 1;
        end
        
        
      
       
        string;
        
        % if there are RMAX*CMAX labels then make pdf
        if mod(e,RMAX*CMAX) == 0
            fprintf(['placing page footer @:' num2str(e) ' page.\n']);
            fprintf(fid,pageFooterString);
            %fprintf(fid,'\\newline\n');
            %{
            toCo{end+1} = ['page' num2str(cnt) '.pdf'];
            print(gcf,'-dpdf',['page' num2str(cnt) '.pdf']);
            cnt = cnt + 1;
            close all
            
            %}
            figSaved = 1;
        end
        fprintf(['end generating QR label latex code:' num2str(e) '\n']);
    end
    
    
    
    
    % save last sheet if it is not full
    if ~figSaved
        %fprintf(fid,pageFooterString);
        %{
        toCo{end+1} = ['page' num2str(cnt) '.pdf'];
        print(gcf,'-dpdf',['page' num2str(cnt) '.pdf']);
        cnt = cnt + 1;
        close all
        figSaved = 1;
        %}
    end
    
    fprintf(fid,documentFooterString);
    fclose(fid);
    CMD = ['pdflatex --interaction=nonstopmode ' oPath 'masterDoc.tex'];
    system(CMD),'-echo';
    drawnow
    close all
end

%{
    csvFile = '/home/nate/Documents/EDGAR.csv';
    csvFile = '/home/nate/Downloads/MBS 2016 QR code sheet-1.csv';
    csvFile = '/home/nate/Downloads/Mike.csv';
    csvFile = '/home/nate/Downloads/Anibus.csv';
    csvFile = '/home/nate/Downloads/density QR.full.csv';
    csvFile = '/home/nate/Documents/EDGAR.csv';
    %csvFile = '/home/nate/Downloads/qr_csv.csv - qr_csv.csv.csv';
    csvFile = '/home/nate/Downloads/2017rep2pop45QRs.csv';
    genQRsheets_latex(csvFile,'./output/');
    csvFile = '/home/nate/Downloads/Trevor.csv';
    genQRsheets_latex(csvFile,'./Trevor_output2/');

%}