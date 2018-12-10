function [doc] = JSONcompile(FileList,oPath)
    % make the output table
    %doc = table;
    cnt = 0;
    mkdir(oPath);
    
   
    
    for e = 1:numel(FileList)
        try
          
            
            tic
            % reporting
            fprintf(['start loading: ' num2str(e) ':' num2str(numel(FileList)) '\n']);
            string = fileread(FileList{1});
            tmpD = jsondecode(string);
            %{
            % load the json file
            tmpD = loadjson(FileList{e});
            %}
            
            
            % get file parts
            [tmpP tmpNM tmpEXT] = fileparts(FileList{e});
            % get the name
            tmpNM = tmpNM(1:end-5);
            % get the fields from the json 
            f1 = fields(tmpD);
            % get the document
            tmpD = tmpD.(f1{1});
            % for each object in the document
            
            
            for object = 1:numel(tmpD)
                % get the fields
                f2 = fields(tmpD(object));
                % increase the object count
                cnt = cnt + 1;
                % for each field
                for f = 1:numel(f2)
                    % get the type
                    type = tmpD(object).(f2{f}).type;
                    
                    if ~strcmp(type,'char')
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % get from file
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % get the size of the data object
                        sz = tmpD(object).(f2{f}).size;
                        % get the data
                        data = tmpD(object).(f2{f}).data;
                        % get the human name
                        name = tmpD(object).(f2{f}).name;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % put to table
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        doc.fileName{cnt,1} = tmpNM;
                        doc.objectNumber{cnt,1} = object;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        if prod(sz) == 1
                            sz = 1;
                        end
                        if ndims(sz) > numel(name)
                            name{2} = name{1};
                            name{1} = 'default';
                        end

                        ind = cell(1,ndims(sz));
                        for d = 1:prod(sz)
                            [ind{:}] = ind2sub(sz,d);
                            header = '';
                            for n = 1:numel(ind)
                                header = [header char(name{n}) '__' num2str(ind{n}) '_'];
                            end
                            header(end) = [];
                            doc.(header){cnt,1} = data(d);
                        end
                    end

                end
                fprintf(['done with :' num2str(object) ':' num2str(numel(tmpD)) '@' num2str(toc) '\n']);
            end
            fprintf(['done with :' num2str(e) ':' num2str(numel(FileList)) '@' num2str(toc) '\n']);
        catch ME
            getReport(ME)
            fprintf(['failed with :' num2str(e) ':' num2str(numel(FileList)) '@' num2str(toc) '\n']);
        end
    end
    doc = struct2table(doc);
    fprintf(['start write table : \n']);tic
    writetable(doc,[oPath filesep date '_JSON_compiledData.csv'])
    fprintf(['stop write table :' num2str(toc) '\n']);
end

%{
    %massDownload('/iplant/home/kmichel/maizeData/return/cobData/', '.json','/home/nate/Download/testCob/')
    FilePath = '/home/nate/Download/testCob/';
    FileList = {};
    FileExt = {'json'};
    FileList = gdig(FilePath,FileList,FileExt,1);

    FilePath = '/home/nate/Download/testCob/';
    FileList = {};
    FileExt = {'json'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    FileList{1} = '/home/nate/Downloads/E_TS_01_jdoc.json';
    FileList{1} = '/home/nate/Downloads/202380_imageData_ear1_jdoc.json';
    FileList{1} = '/home/nate/output/{Plot_874}{Experiment_42}{Planted_4-18-2016}{SeedSource_BAM349-3}{SeedYear_2016}{Genotype_B73}{Treatment_4x 0C-6C stress}{PictureDay_15}_jdoc.json';
    doc = JSONcompile(FileList,'./output/')


    FilePath = '/home/nate/Downloads/caliOSG/';
    FileList = {};
    FileExt = {'json'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    doc = JSONcompile(FileList,'/home/nate/Downloads/caliOSG/');


    FilePath = '/home/nate/Downloads/taraJSON/';
    FileList = {};
    FileExt = {'json'};
    FileList = gdig(FilePath,FileList,FileExt,1);
    doc = JSONcompile(FileList,'/home/nate/Downloads/JUNKY/');


    R = readtext('/home/nate/Downloads/4-18-18cobJSON2016');
    FileList = R(2:5);
    doc = JSONcompile(FileList,'/home/nate/Downloads/JUNKY/');

    R = readtext('/home/nate/Downloads/4-18-18cobJSON2016');
    R = R(1:5);


    jfile = {'/home/nate/Downloads/6_jdoc.json'};
    doc = JSONcompile(jfile,'/home/nate/Downloads/JUNKY/');






%}