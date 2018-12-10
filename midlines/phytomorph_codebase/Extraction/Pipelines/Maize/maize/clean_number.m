FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies/';
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% 1) sort SET
for e = 1:numel(SET)
    numel(SET{e});
    if numel(SET{e}) > 61
        fn = {SET{e}{62:end}};
        for i = 1:numel(fn)
            f = fn{i}
            delete(f);
        end
    end
end
%% rename folders
convertFile = '/home/nate/Downloads/convertShort.csv';
[convertList]= readtext(convertFile);
import java.util.HashMap;
O2N = HashMap();
for e = 1:size(convertList,1)
    O2N.put(convertList{e,1},convertList{e,2});
end

%%
import java.io.File
java.io.File('/home/nate/Downloads/convertShort.csv').renameTo(java.io.File('/home/nate/Downloads/convertShort2.csv'))
%%
java.io.File('/home/nate/Downloads/convertShort2.csv').renameTo(java.io.File('/home/nate/Downloads/convertShort.csv'))
%% get sample
O2N.get('S5_2067B_E8_E4_E1_D7_D5_D1')
O2N.get('S1_2050B_B7_B6_B5_B2_A4')
%% loop over old
missing = 0;
for e = 1:numel(SET)
    [p n ex] = fileparts(SET{e}{1});
    fidx = strfind(p,filesep);
    k = p(fidx(end)+1:end);
    v = O2N.get(k);
    if isempty(v)
        missing = missing + 1;
    else
        NEWname = [p(1:fidx(end)-1) filesep v];
        OLDname = p;
        fprintf([k '-->' v '\n']);
        fprintf('*****\n');
        fprintf([OLDname '\n']);
        fprintf([NEWname '\n']);
        fprintf('*************************\n');
        java.io.File(OLDname).renameTo(java.io.File(NEWname));
        %java.io.File(NEWname).renameTo(java.io.File(OLDname));
    end
end