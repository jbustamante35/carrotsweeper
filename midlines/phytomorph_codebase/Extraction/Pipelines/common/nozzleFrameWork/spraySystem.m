function [] = spraySystem(FileList,nozzleName)

     SUBimageSource = dataSource(FileList,@(X)double(imread(X)),1);
end