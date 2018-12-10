function [] = genQRLargeFormatSheets(txtFile,oPath)
    mkdir(oPath);
    CMD = ['/de-app-work/generate_qrcode_R_script.pl --in ' txtFile ' --out ' oPath 'QRcode.R'];
    [status,cmdout] = system(CMD,'-echo');
    CMD = ['R -f ' oPath 'QRcode.R'];
    [status,cmdout] = system(CMD,'-echo');
end