function [I] = generalizeLoader(fileName,loaderType,loaderArgs)
    switch loaderType
        case 'stomata_histonormalize'
            I = imread(fileName);
            I = imhistmatch(I,loaderArgs{1});
    end
end