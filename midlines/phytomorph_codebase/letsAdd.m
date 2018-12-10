function [c,m] = letsAdd(a,b)
    c = a + b;
    m = a - b;
end

%{

    result = letsAdd(1,2);

    cnt = 1;
    for a = 1:10
        for b = 1:10
            r(cnt) = letsAdd(a,b);
            cnt = cnt +1;
        end
    end


    func = cFlow('letsAdd');
    func.setMCRversion('v920');
    cnt = 1;
    for a = 1:2
        for b = 1:2
            [c{cnt} m{cnt}] = func(a,b);
            cnt = cnt + 1;
        end
    end


    auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
    auth = auth{1};
    func.submitDag(auth,50,50);





    [o1] = cFlowLoader(c{1});
    I = imread('/home/nate/Downloads/GH3_760_Oh43_100.tif');
    I = imread('/iplant/home/mrwhite4/maizeData/earData/Scan3-170213-0072.tif');


%}