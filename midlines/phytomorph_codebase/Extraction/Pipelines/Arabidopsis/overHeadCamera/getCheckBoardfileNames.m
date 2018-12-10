function [fileNames] = getCheckBoardfileNames(table)
    fidx = table.checkerBoard == true;
    fileNames = table.fileName(fidx);
end