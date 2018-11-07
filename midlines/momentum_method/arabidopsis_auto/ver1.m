%% find files and generate rand number
pth = '/mnt/scratch3/users/nmiller/phytoMorph/outPort/'
dataFiles = gdig(pth,{},{'mat'},1);
numel(dataFiles);

%% make radial patch
[theta,rb] = generate_FFT_patch(11);


%%
for e = 1%:10
    load(dataFiles{e})
    rootN = 1;
    I = imread(d{rootN}{1}.image.fileName);        
    % handle flip
    [I direc] = handleFLIP(I,d{rootN}{1}.flipD);
    % resize I
    I = imresize(I,.25);
    % preallocate results
    f = zeros([size(I,1) size(I,2) size(rb,3)]);
    % conv
    I = double(I);
    for i = 1:size(rb,3)
        f(:,:,i) = conv2(I,rb(:,:,i),'same');
    end
    % reshape
    fsz = size(f);
    f = reshape(f,[size(f,1)*size(f,2) size(f,3)]);

    % get magnitude
    f = abs(f);
    h = ones(1,10);
    
    
    % loop over each feature vector
    df = imdilate(f,h);
    for i = 1:size(f,1)
        % nonmax supression
        %df = imdilate(f(i,:),h);
        fidx = find(f(i,:) == df(i,:));
        dff = f(i,fidx);
        [dff sidx] = sort(dff,'descend');
        % get the first N=2 frequecies and magntudes
        fp(:,:,i) = [fidx(sidx(1:2));dff(1:2)];
        i
    end
    
    fp = reshape(fp,[4 size(fp,3)]);
    K = reshape(fp(4,:),fsz(1:2));
    imshow(K,[]);
    
    [J sidx] = sort(abs(f),3,'descend');
    
    
    
    for i = 1:size(rb,3)
        f(:,:,i) = f(:,:,i) * max(max(abs(f(:,:,i))))^-1;
    end
    
    for i = 1:size(rb,3)
        imshow(abs(f(:,:,i)),[]);
        drawnow
    end
end
%%
Tr.learn();
%%
Tr.viewManifold();
%% patch

h = figure;
patchS.view(h,[],[]);

%% look at curves

close all

Tr.viewCoandDo();
close all




