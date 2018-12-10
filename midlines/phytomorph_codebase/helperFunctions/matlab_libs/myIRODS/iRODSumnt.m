function [] = iRODSumnt(mntPoint)
    cmd = ['fusermount -u ' mntPoint];
    system(cmd);
end