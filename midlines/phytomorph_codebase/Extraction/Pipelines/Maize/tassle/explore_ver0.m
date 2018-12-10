FilePath = '/mnt/spaldingdata/nate/tips/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% loading
for e = 1:numel(FileList)
    tmp = load(FileList{e});
    DM(:,:,e) = tmp.spikeTip;
    e
end
%% chop DM
DM(end-(100-1):end,:,:) = [];
%%
nDM = bsxfun(@minus,DM,mean(DM,1));
fDM = fft(nDM,[],1);
ang = angle(fDM);
for e = 1:size(ang,2)
    for e2 = 1:size(ang,3)
        if ang(1,e,e2) <= 0 
            ang(:,e,e2) = ang(:,e,e2) + 2*pi; 
        end
    end
end
ang = unwrap(ang,[],1);
amp = abs(fDM);
X = real(fDM);
Y = imag(fDM);
%%
close all
plot(X(:,100,1),Y(:,100,1));
figure;
plot(X(:,100,1))
%%
syn = amp.*exp(i*ang);
syn1 = ifft(syn,[],1);
close all
imshow(syn1(:,:,1),[]);
%% pca along width only
for w = 1:size(amp,2)
    [aS(:,w,:) aC(:,w,:) aU(:,w) aE(:,:,w)] = PCA_FIT_FULL_T(squeeze(amp(:,w,:)),3);
end
%% perfect recon
perfSyn = amp.*exp(i*ang);
perfImg = ifft(perfSyn,[],1);

perfImg2 = ifft(fDM,[],1);
%% syn recon
close all
syn = aS.*exp(i*ang);
syn1 = ifft(syn,[],1);
syn1 = bsxfun(@plus,syn1,mean(DM,1));
imshow(syn1(:,:,10),[]);
%%

%% pca along width only
clear xxS xxC xxU xxE yyS yyC yyU yyE
for w = 1:size(X,2)
    [xxS(:,w,:) xxC(:,w,:) xxU(:,w) xxE(:,:,w)] = PCA_FIT_FULL_T(squeeze(X(:,w,:)),10);
    [yyS(:,w,:) yyC(:,w,:) yyU(:,w) yyE(:,:,w)] = PCA_FIT_FULL_T(squeeze(Y(:,w,:)),10);
end
%% simulate
syn = xxS + i*yyS;
%syn = X + i*Y;
syn1 = ifft(syn,[],1);
%syn1 = bsxfun(@plus,syn1,mean(DM,1));
close all
imshow(syn1(:,:,10),[]);

%% try whole decompose of amp
szAMP = size(amp);
vec_amp = reshape(amp,[szAMP(1)*szAMP(2) szAMP(3)]);

[vS vC vU vE] = PCA_FIT_FULL_T(vec_amp,3);

szANG = size(ang);
vec_ang = reshape(ang,[szANG(1)*szANG(2) szANG(3)]);

[aS aC aU aE] = PCA_FIT_FULL_T(vec_ang,5);
vS = reshape(vS,szAMP);
aS = reshape(aS,szANG);
syn = vS.*exp(i*ang);
%syn = vS.*exp(i*aS);
%syn = amp.*exp(i*ang);
syn1 = ifft(syn,[],1);
%%
szX = size(X);
vec_X = reshape(X,[szX(1)*szX(2) szX(3)]);
[xS xC xU xE] = PCA_FIT_FULL_T(vec_X,3);
szY = size(Y);
vec_Y = reshape(Y,[szY(1)*szY(2) szY(3)]);
[yS yC yU yE] = PCA_FIT_FULL_T(vec_Y,3);
xS = reshape(xS,szX);
yS = reshape(yS,szY);
%%
close all
plot3(xC(1,:),xC(2,:),xC(3,:),'.');
figure
plot3(yC(1,:),yC(2,:),yC(3,:),'r.');
%% simulate
syn = xS + i*yS;
%syn = X + i*Y;
syn1 = ifft(syn,[],1);
syn1 = bsxfun(@plus,syn1,mean(DM,1));
%%
close all
plot3(vC(1,:),vC(2,:),vC(3,:),'.')
%% try to correct phase
for e = 1:size(ang,3)
    for w = 1:size(ang,2)
        for r = 1:size(ang,1)-1
            
            cur = ang(r,w,e);
            nxt = ang(r+1,w,e);
            %{
            if nxto < 0
                nxto = pi + (pi - abs(nxt));
            else
                nxto = -pi - (pi - nxt);
            end
            choi = [(-10:1:10)*2*pi + nxt (-10:1:10)*2*pi + nxto ];
            %choi = [nxt nxto];
            %}
            choi = [-100:1:100]*2*pi + nxt;
            [v sidx] = min(abs(cur - choi));
            ang(r+1,w,e) = choi(sidx);
        end
    end
    e
end
%% align vecs
template = ang(:,30,1);
for e = 1:size(ang,2)
    for e2 = 1:size(ang,3)
        if ang(:,e,e2)'*template < 0
            ang(:,e,e2) = -ang(:,e,e2);
        end
    end
    e
end