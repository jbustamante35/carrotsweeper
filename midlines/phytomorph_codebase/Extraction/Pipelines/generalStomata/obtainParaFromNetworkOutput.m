function [para] = obtainParaFromNetworkOutput(networkOutput)
    para = [networkOutput(9:11) 0 0];
end