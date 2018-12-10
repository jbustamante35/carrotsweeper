function [STORE,DIS] = findProgram(disp,ML,LOOP,Y,X,MUTATE_NUMBER,MAX_WID,BYTE_SIZE,PROGRAM)

    if disp
        h1 = figure;
        h2 = figure;
    end
    
    eP = [];
    delta = [];
    for ml = 1:ML
        
        %{
        vars = 4^2;
        degreeFreedom = nchoosek(16,0);
        degreeFreedom = nchoosek(8,1);
        degreeFreedom = nchoosek(4,2);
        degreeFreedom = nchoosek(2,3);
        degreeFreedom = nchoosek(1,4);
        degreeFreedom = nchoosek(0,5);
        %}
        PROGRAM_TEST = randi(255,1,MAX_WID,'uint8');

        for loop = 1:LOOP
            if nargin > 8
                delta(ml,loop) = bitwise_hamming(PROGRAM_TEST',PROGRAM');
                figure(h1);
                plot(delta');
            else
                delta(ml,loop) = 0;
            end
            %delta(ml,loop) = norm(double(PROGRAM_TEST) - double(PROGRAM));



            toMut = PROGRAM_TEST;
            cnt = 1;
           

            NUM_BYTE_TO_MUT = min(MUTATE_NUMBER,size(PROGRAM_TEST,2));
            rand_MUT_byte = randi(size(toMut,2),NUM_BYTE_TO_MUT,1);
            progB = zeros(NUM_BYTE_TO_MUT*(2^BYTE_SIZE),size(PROGRAM_TEST,2),'uint8');
            
            %{
            %fprintf(['Starting mutation of programs @ bytes @ ' num2str(NUM_BYTE_TO_MUT) '.\n']);
            tm = clock;
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
            ftm = etime(clock,tm);
            %fprintf(['Starting mutation of programs @ bytes @  ' num2str(ftm/60) ' min.\n']);
            %}
            
            
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
            



            
            %fprintf(['Starting eval of programs @ ' num2str(size(progB,1)) '.\n']);
            tm = clock;
            preY = bo(X,progB(1,:));
            ftm = etime(clock,tm);
            %fprintf(['Estimated delta-T @ ' num2str(size(progB,1)*ftm/60) ' min.\n']);
            tm = clock;


            
            preY = boo(X,progB);
            ME = zeros(size(progB,1),1);
            for iter = 1:size(progB,1)
                ME(iter) = -mutualinfo(preY(:,iter),Y);
                %ME(iter) = bitwise_hamming(preY(:,iter),Y);
            end
            
            

            %{
            ME = zeros(size(progB,1),1);
            for iter = 1:size(progB,1)
                preY = bo(X,progB(iter,:));
                if ~all(preY == preY2(:,iter))
                    break
                end
                ME(iter) = mutualinfo(uint8(preY),Y);
            end
            %}

            ftm = etime(clock,tm);
            %fprintf(['Starting eval of programs @' num2str(size(progB,1)) '.\n']);
            %fprintf(['Real delta-T @ ' num2str(ftm/60) ' min.\n']);
            

            initY = bo(X,toMut);
            %initE = mutualinfo(uint8(initY),Y);


            %[finalE,mi] = max(ME);
            [finalE,mi] = min(ME);
            
            PROGRAM_TEST = progB(mi,:);
            
            
            finalE = bitwise_hamming(boo(X,PROGRAM_TEST),Y);
            eP(ml,loop) = finalE;

            if disp
                figure(h2);
                plot(eP');
                drawnow
            end
        end
        STORE(ml,:) = PROGRAM_TEST;
        DIS(ml,1) = eP(ml,end);
        DIS(ml,2) = delta(ml,end);
    end
    if disp
        close(h1)
        close(h2)
    end
end