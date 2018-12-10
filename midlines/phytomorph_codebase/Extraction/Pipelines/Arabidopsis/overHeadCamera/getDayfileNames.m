function [fileNames] = getDayfileNames(table)
    fidx = table.day == true & table.checkerBoard == false;
    fileNames = table.fileName(fidx);
end