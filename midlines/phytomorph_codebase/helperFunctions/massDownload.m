function [] = massDownload(sourceLocation,searchString,targetLocation)
    mkdir(targetLocation)
    %CMD = ['ils ' sourceLocation ' | grep ' searchString ' | xargs -L 1 -P 20 sh -c ''iget -f "' sourceLocation '"${0}" ''' targetLocation ''];
    %CMD = ['ils ' sourceLocation ' | grep ' searchString ' | sed -e "s/^*//" | xargs -d ''\n'' -L 1 -P 20 sh -c ''iget -f "' sourceLocation '"${0}" ' targetLocation ''''];
    %CMD = ['ils ' sourceLocation ' | grep ' searchString ' | sed -e "s/^[ \t]*//" | xargs -L 1 -P 10 sh -c ''iget -f "' sourceLocation '${0}" "' targetLocation  '"'''];
    
    %CMD = ['ils ' sourceLocation ' | grep ' searchString ' | sed -e "s/^[ \t]*//"'];
    
    CMD = ['ils ' sourceLocation ' | grep ' searchString ' | xargs -d ''\n'' -L 1 -P 20 sh -c ''sed -e "s/^[ \t]*//" | iget -f "' sourceLocation '"${0}" ' targetLocation ''''];
    
    %ils | grep json | sed -e "s/^[ \t]*//" | xargs -L1 sh -c 'iget -f ${0} /home/nate/'
    [status,cmdout] = system(CMD,'-echo');
end
%{
    massDownload('/iplant/home/kmichel/maizeData/return/cobData/', '.json','/home/nate/Download/testI/');
    massDownload('/iplant/home/petersoapes/returns/return5/', '.mat','/home/nate/Downloads/learnDNA/');
    

    massDownload('/iplant/home/canibas/fastOSG/', '.mat','/mnt/tetra/nate/caliSample/');
    list = readtext('/home/nate/Downloads/All_030218');

    massDownload('/iplant/home/canibas/fastOSG/', '.tif','/mnt/tetra/nate/caliSampleRAW/');
    massDownload('/iplant/home/canibas/maizeData/seedlingData/', '.nef','/mnt/tetra/nate/caliSampleNEF/');
    list = readtext('/home/nate/Downloads/All_030218');


    massDownload('/iplant/home/phytomorphuser/workITOUT_FAST/', '.csv','/mnt/tetra/nate/forLeakey/');

    dataPath = ['/iplant/home/hirsc213/maizeData%'];

    massDownload('/iplant/home/hirsc213/maizeData/', '.nef','/mnt/tetra/nate/seedlingDATApile/');

    massDownload('/iplant/home/ssu1/2017GWAS/Output/', '.mat','/mnt/tetra/nate/MATforshiheng/');







%}
