function [z] = bo(a,b)
    %tic
    z = any(bitand(a,b) > 0,2);
    %toc
    %{
    tic
    z = ~all(bitand(a,b) == 0,2);
    toc
    %}
    %{
    tic
    z = sum(bitand(a,b),2)>0;
    toc
    %}
    %z = uint8(0);
    %{
    fidx = [find(c)];
    for e = 1:numel(fidx)
        z = bitset(z,fidx(e));
    end
    %}
end


%{
    % simple
    fidx = [];
    while isempty(fidx)
        N = 1400;
        W = 80;
        GN = 2;
        A = uint8(randi(255,N,W));
        gate = uint8(randi(255,GN,W));

        tmi = clock;
        Aoo = boo(A,gate);
        tm = etime(clock,tmi);
        fprintf(['Done BOO:' num2str(tm) '\n']);

        tmi = clock;
        Ao = zeros(size(A,1),size(gate,1));
        for g = 1:size(gate,1)
            Ao(:,g) = bo(A,gate(g,:));
        end
        tm = etime(clock,tmi);
        fprintf(['Done BO:' num2str(tm) '\n']);


        fidx = find(double(Aoo(:))~=Ao(:));
    end


    A = uint8(randi(255,N,W));
    gate = uint8(randi(255,1,W));





   
    BYTE_SIZE = 3;                              % log size of a byte 
    MAX_B = 8;                                  % max bits to handle at once
    MAX_B = 12;                                 % max bits to handle at once
    MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide



    TALL = 10000;


    D = zeros(TALL,MAX_WID,'uint8');

    INT = double(randi((2^MAX_B)-1,TALL,1));    % encode this- int to test no zeros!
    BYTE = floor(double(INT)*(2^BYTE_SIZE)^-1)+1;
    REM = rem(INT-1,2^BYTE_SIZE)+1;
    POS = sub2ind(size(D),(1:size(D,1))',BYTE);
    D(POS) = bitset(D(POS),REM);





    % sparse speed test
    sD = sparse((1:size(D,1))',double(INT),1,size(D,1),MAX_WID*2^BYTE_SIZE);
    




    %% test idea
    xD = zeros(TALL,MAX_WID,'uint8'); % pre allocate data
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
    [STORE,DIS] = findProgram(true,50,50,Y,X,100,MAX_WID,BYTE_SIZE,PROGRAM);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% TEST SQUEEZE BITS
    lemon = randi(2,100,3,'uint8');
    [juice] = squeezeQbits(lemon,2,5,'uint32');



    [Y] = unpackQbits(X,2);
    cY = Y'*Y;



    % sparse speed test for program
    cnt = 1;
    for byte = 1:size(PROGRAM,2)
        for bit = 1:2^BYTE_SIZE
            sP(cnt) = bitget(PROGRAM(byte),bit);
            cnt = cnt + 1;
        end
    end
    sPROGRAM = sparse(1:numel(sP),1,double(sP));


    %{
    % random program tries
    tic
    for iter = 1:10000000
        PROGRAM_TEST = randi(255,1,MAX_WID,'uint8');
        preY = bo(X,PROGRAM_TEST);
        ME(iter) = mutualinfo(uint8(preY),Y);
    end
    toc
    %}
    


    close all
    h1 = figure;
    h2 = figure;
    eP = [];
    delta = [];
    for ml = 1:10
        PROGRAM_TEST = randi(255,1,MAX_WID,'uint8');

        for loop = 1:2000
            delta(ml,loop) = norm(double(PROGRAM_TEST) - double(PROGRAM));

            figure(h1);
            plot(delta');


            toMut = PROGRAM_TEST;
            cnt = 1;
           

            NUM_BYTE_TO_MUT = min(800,size(PROGRAM_TEST,2));
            rand_MUT_byte = randi(size(toMut,2),NUM_BYTE_TO_MUT,1);
            progB = zeros(NUM_BYTE_TO_MUT*(2^BYTE_SIZE),size(PROGRAM_TEST,2),'uint8');
            %{
            % for stochastic gradient decent - mut each bit in rand byte - type I
            for b = 1:NUM_BYTE_TO_MUT
                byte = rand_MUT_byte(b);
                for bit = 1:(2^BYTE_SIZE)
                    newP = toMut;
                    newP(byte) = bitset(newP(byte),bit,not(bitget(newP(byte),bit)));
                    progB(cnt,:) = newP;
                    cnt = cnt + 1;
                end
            end
            %}
            
            %{
            % for gradient decent
            for byte = 1:MAX_WID
                for bit = 1:(2^BYTE_SIZE)
                    newP = toMut;
                    newP(byte) = bitset(newP(byte),bit,not(bitget(newP(byte),bit)));
                    progB(cnt,:) = newP;
                    cnt = cnt + 1;
                end
                byte;
                MAX_WID;
            end
            %}



            
            fprintf(['Starting eval of programs @ ' num2str(size(progB,1)) '.\n']);
            tm = clock;
            preY = boo(X,progB(1,:));
            ftm = etime(clock,tm);
            fprintf(['Estimated delta-T @ ' num2str(size(progB,1)*ftm/60) ' min.\n']);
            tm = clock;


            
            preY = boo(X,progB);
            ME = zeros(size(progB,1),1);
            for iter = 1:size(progB,1)
                ME(iter) = mutualinfo(preY(:,iter),Y);
            end
            

            %{
            ME = zeros(size(progB,1),1);
            for iter = 1:size(progB,1)
                preY = bo(X,progB(iter,:));
                ME(iter) = mutualinfo(uint8(preY),Y);
            end
            %}

            ftm = etime(clock,tm);
            fprintf(['Starting eval of programs @' num2str(size(progB,1)) '.\n']);
            fprintf(['Real delta-T @ ' num2str(ftm/60) ' min.\n']);
            

            initY = bo(X,toMut);
            initE = mutualinfo(uint8(initY),Y);


            [finalE,mi] = max(ME);
            PROGRAM_TEST = progB(mi,:);
            eP(ml,loop) = finalE;


            figure(h2);
            plot(eP');
            drawnow
        end
    end


    




    

    %{
    for e = 1:TALL
        D(e,BYTE(e)) = bitset(D(e,BYTE(e)),REM(e));
    end
    %}


%}