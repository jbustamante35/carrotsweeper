function [] = viewer(fn)
    list = hdf5ls(fn);
    cnt = 1;
    
    dataStore = {};
    for g = 1:numel(list.Groups)
        dataStore
        if ~isempty(list.Groups(g).meta)
            
            flds = fields(list.Groups(g).meta);
            
            
            %{
            if any(strcmp(flds,'type'))
                
                if strcmp(list.Groups(g).meta.type,'image');
                    
                    external_uriList.{cnt} = list.Groups(g).meta.uri;
                    internal_uriList{cnt} = list.Groups(g).uri;
                    cnt = cnt + 1;
                    
                end
                
                
                if strcmp(list.Groups(g).meta.type,'curve');
                    
                    external_uriList{cnt} = list.Groups(g).meta.uri;
                    internal_uriList{cnt} = list.Groups(g).uri;
                    cnt = cnt + 1;
                    
                end
                
            end
            %}
        end
    end
    
    %%% sort list
    toSort = 1;
    if toSort
        for fl = 1:numel(external_uriList)
            [pth,nm{fl},ext] = fileparts(external_uriList{fl});
            num_nm(fl) = str2num(nm{fl});            
        end
        [junk sidx] = sort(num_nm);
        external_uriList = external_uriList(sidx);
        internal_uriList = internal_uriList(sidx);
    end
    
    %%% get tensors for attached to images    
    for i = 1:numel(external_uriList)
        
    end

    %%% view
    for fl = 1:numel(uriList)        
        I = imread(uriList{fl});
        imshow(I,[])
        drawnow
    end
    
    %{

    % for each contour and midline - there may be many
    for obj = 1:numel(dataStore.img{fl}.objects{obj})
        for nten = 1:numel(pidx)
            plot(dataStore.img{fl}.objects{obj}.tensors{pidx(nten)}.data(:,2),dataStore.img{fl}.objects{obj}.tensors{pidx(nten)}.data(:,1));
            LEG{nten} = dataStore.img{fl}.objects{obj}.tensors{pidx(nten)}.name;
        end
    end           
            
            % disp the quiver
            pidx0 = [8];
            pidx1 = [10];
            % for each contour and midline - there may be many
            for obj = 1:numel(dataStore.img{fl}.objects{obj})
                for nten = 1:numel(pidx0)                                        
                    quiver(dataStore.img{fl}.objects{obj}.tensors{pidx0(nten)}.data(:,2),...
                           dataStore.img{fl}.objects{obj}.tensors{pidx0(nten)}.data(:,1),...
                           dataStore.img{fl}.objects{obj}.tensors{pidx1(nten)}.data(:,2,1),...
                           dataStore.img{fl}.objects{obj}.tensors{pidx1(nten)}.data(:,1,1));
                       
                   quiver(dataStore.img{fl}.objects{obj}.tensors{pidx0(nten)}.data(:,2),...
                       dataStore.img{fl}.objects{obj}.tensors{pidx0(nten)}.data(:,1),...
                       dataStore.img{fl}.objects{obj}.tensors{pidx1(nten)}.data(:,2,2),...
                       dataStore.img{fl}.objects{obj}.tensors{pidx1(nten)}.data(:,1,2));
                    LEG{nten} = dataStore.img{fl}.objects{obj}.tensors{pidx1(nten)}.name;
                end
            end
            
            
            legend(LEG);
            
            
            
        end
    %}
    end
  