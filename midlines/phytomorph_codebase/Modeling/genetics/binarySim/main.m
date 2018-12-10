%% test idea
BYTE_SIZE = 3;                              % log size of a byte 
MAX_B = 8;                                  % max bits to handle at once
%MAX_B = 12;                                 % max bits to handle at once
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide
TALL = 10000;
xD = zeros(TALL,MAX_WID,'uint8'); % pre allocate data,'uint8'); % pre allocate data
%INT = uint8(randi((2^MAX_B)-1,TALL,1)); % int to test no zeros!! % generate random data
INT = double(randi((2^MAX_B)-1,TALL,1)); % int to test no zeros!! % generate random data
BYTE = floor(double(INT)*(2^BYTE_SIZE)^-1)+1;   
REM = rem(INT-1,2^BYTE_SIZE)+1;
POS = sub2ind(size(D),(1:size(D,1))',BYTE);
xD(POS) = bitset(D(POS),REM);

X = xD;
bit1 = randi(8,1);
bit2 = randi(8,1);
PROGRAM = randi(255,1,MAX_WID,'uint8');

%PROGRAM(1:end-1) = 0;

Y = uint8(bo(X,PROGRAM));

[STORE,DIS] = findProgram(10,2000,Y,X,10,MAX_WID,BYTE_SIZE,PROGRAM);


%%
mex -largeArrayDims -lpthread boo.cpp
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TALL = 100000;    
MAX_B = 8;
BYTE_SIZE = 3;                                                          % log size of a byte
MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);
TALL = 1000;
[INT] = samplePureStates(MAX_B,TALL);


MAX_B = 4;
[INT] = samplePureStates(MAX_B,TALL);
[INT] = squeezeQbits([INT INT],4,3,'uint8');
MAX_B = 8;


[X] = quantumEncode(MAX_B,INT);
PROGRAM = randi(255,1,MAX_WID,'uint8');
Y = boo(X,PROGRAM);

iterations = 200;
repeats = 10;
[STORE,DIS] = findProgram(true,repeats,iterations,Y,X,100,MAX_WID,BYTE_SIZE,PROGRAM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
