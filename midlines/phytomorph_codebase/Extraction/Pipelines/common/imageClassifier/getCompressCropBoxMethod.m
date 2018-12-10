function [U,E] = getCompressCropBoxMethod(toCompress,numC)
    % stack and compress
    XStack = zeros([size(toCompress{1},1) size(toCompress{1},2)*numel(toCompress)]);
    str = 1;
    jmp = size(toCompress{1},2);
    stp = str + jmp - 1;
    for e = 1:numel(toCompress)
        XStack(:,str:stp) = toCompress{e};
        str = stp + 1;
        stp = str + jmp - 1;
    end
    [~,~,U,E] = PCA_FIT_FULL_T(XStack,numC);
end