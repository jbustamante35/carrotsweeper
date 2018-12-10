function [] = surKur_loop(patchStack,para)
    % number of trial
    num_tr = size(patchStack,ndims(patchStack));
    
    % patch size
    patchSize = size(patchStack);
    patchSize = patchSize(1:end-1);
    % loop over
    for e = 1:num_tr 
        [sK kVec] = surKur(patchStack(:,:,e),para);
    end
end

%{


imshow(patchStack(:,:,e),[]);
imshow(eye(200),[]);
hold on
[n1 n2] = ndgrid(1:1:size(patchStack,1),1:1:size(patchStack,2));
[n1 n2] = ndgrid(1:1:200,1:1:200);

quiver(n2,n1,kVec(:,:,1),kVec(:,:,2),'b');
quiver(n2,n1,kVec(:,:,3),kVec(:,:,4),'r');

quiver(n2,n1,kVec(:,:,2),kVec(:,:,1),'g');
quiver(n2,n1,kVec(:,:,4),kVec(:,:,3),'m');



[x y V] = impixel(patchStack(:,:,e)/255);
[P] = flow(kVec(:,:,3:4),[x y],500,.2,patchStack(:,:,e));

%}