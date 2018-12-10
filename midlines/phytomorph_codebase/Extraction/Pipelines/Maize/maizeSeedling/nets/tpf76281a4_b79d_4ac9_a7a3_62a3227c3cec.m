function [Y,Xf,Af] = tpf76281a4_b79d_4ac9_a7a3_62a3227c3cec(X,~,~)
%TPF76281A4_B79D_4AC9_A7A3_62A3227C3CEC neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 29-Aug-2017 22:01:18.
% 
% [Y] = tpf76281a4_b79d_4ac9_a7a3_62a3227c3cec(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timsteps
%   Each X{1,ts} = 20xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

  % ===== NEURAL NETWORK CONSTANTS =====
  
  % Input 1
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020144;-780.471896742374;-460.20555521635;-234.31060487846;-540.291178237133;-572.142715533789;-192.819908210792;-263.800490630852;-210.820779921779;-258.183122033521;-286.590370604323;-169.741999038737;-271.750472270932;-218.52158935392;-373.104709120361;-120.218367860414;-260.914631799835;-176.115238468635];
  x1_step1_gain = [0.000342540442445779;0.00027748074658588;0.000806013241167019;0.00156677380853905;0.00168393926281428;0.00269046527889069;0.00198084992045411;0.00236361816544284;0.00245784730269611;0.00195235680232596;0.00282106579415896;0.00424658079482962;0.0031726176926076;0.005554586874619;0.00429724736478088;0.00484270636687461;0.00333931608099047;0.00598913897545103;0.00384783002604039;0.00642488839732223];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-0.7436055247597820328;1.5687195182945041783;-1.0940073373418544023;1.3742892976976652797;-0.10955951157601855062;-0.13760961550907935269;2.9442593046474923035;4.110929582319332809;1.9274024442609878349;0.17363483499117482456];
  IW1_1 = [0.25224250241094359204 -0.68889318526338105819 0.082164789731829779296 0.18330261707236703006 -0.09723731576249339259 0.057481450820726552597 -0.19728729868079941667 -0.54825926734564356391 0.78019409138424433703 0.14497745748399704091 0.061412608084417098786 1.0161745125810606183 -0.45814612090477180306 0.52084333725741482901 -0.57688246946294374773 0.3847011577330362031 0.29132039987402325032 -0.07672648363339228117 0.83050542951281458581 0.50325208417846700737;-0.21926633629695005401 -0.30085907893294117077 0.16236793824144693521 0.55384058166443628224 0.36415082120925346132 -0.17787553813458226171 0.1450169728104136746 0.26070894120190818422 -0.64938072954348069477 0.095556288470600256324 -0.34022453797324159197 0.17368123003596991349 0.38468830434597944823 -0.38473112014297300432 -0.32906925812732451542 -0.23366391435797825848 0.58655434811244111515 -0.24611102742756690853 -0.021833582465594011945 0.63072615092877637455;-0.90693868983250613969 0.05511063238819302601 -1.821761848397314365 0.85036453351407936818 -0.19652248535946989505 2.4519125292632599766 -0.37355895262687954794 -0.88061019220303116128 -1.6410674174961508154 2.3930592311435670005 0.89998908705546121212 0.26898620390624722321 0.60414385696556471839 1.0256566861396900769 0.013422677470864131585 -0.076849401799584785011 0.95446914249715575362 0.16141361246900404192 -0.3881045181519735543 -0.5760765537941483494;-1.9018346879083729561 -4.7751548759121265064 -1.5080617727891745261 0.46282969114451910642 -1.118391331672608846 2.3452642058142516568 2.030707686299541237 -0.39015910140297649367 -1.1522065126889720066 -1.2153773372026253075 -2.7332631447407629288 0.45175894819838346983 0.037504396139331375171 0.99532905078650624198 -0.69047633926437979035 -0.92402957155243470311 -4.0215911769320982216 0.054002482641817235731 -1.4716953319142662693 -0.40345373536680934068;0.32766204529821063574 -3.6688227778484017172 2.6127432921164346347 -0.69245677223796797684 0.76733705089814230327 -1.3284542335583862105 0.88292458384404637339 -2.0648087386720437308 -0.1245610620706026489 1.3604539746952866164 -0.35930922296933459625 0.78431670146157961288 -0.30787130900320508475 -2.2299716396013300468 -0.93178420806138106958 0.23471263678839768274 2.4543062486931708222 -1.9031558606916632037 1.1836645340939551119 -0.076424171689879638603;0.86650286462034897816 -2.2101905806928119347 0.64354523436633936928 -1.3763835989967376783 -1.2019990462089491512 1.0315178252172483297 -0.57660690823822857709 1.3484219258811418474 0.51563037257770449884 2.2595684256177994698 -2.1622484521730522822 0.5023015039366975909 -0.11707897375548503849 1.4723617612866781723 1.783346889205093877 -0.65093521408904875258 -1.8532986347231408253 0.049146726006472393755 0.91918592355550154949 1.5742146702388191137;4.768498049279751072 3.4662406253247550758 0.6516004948998113111 -0.0090413517814710644971 2.4370215548010065731 0.79993795645080267853 0.51211493056760104547 -3.0742530384384658682 1.3366055191092904852 -0.75610883572986475887 1.511076379891838517 -0.84134201695257537779 -0.75973341035722197212 0.36341977657540830426 1.2983704464687833013 0.7593037409635557955 -0.71134019849369956034 -1.7689385747634933654 -1.382081367365673108 -1.0917990044211205447;0.30160400579228024975 -0.12107308220896983997 0.061628453579414634034 -0.62304637911516480564 0.98974488199915344744 1.8055519426881274647 1.969600111762311867 -2.4940090119261411061 2.1707789313445307045 3.1674577140548709764 -0.50457328075075358331 1.3381166182499362183 -1.8283814487901484114 -0.394571303735419876 0.52484653848445195212 -0.18613726057206012032 0.81332364681333191303 1.3462052706467124708 2.5882078528720944277 0.10014987138992299809;0.242789124835848108 0.44981276367980349562 0.22559381224972113733 0.12174008013138985784 0.28050617307535530287 0.25559682879054052007 -0.15960933095100018897 0.46916742328371491277 -0.36882506882647736468 -0.55185942742718252507 -0.16386262441614632412 -0.12260578586598237794 -0.58549184886726624022 0.019867441440682045983 0.26139471778310569006 0.42183654986271418297 0.2819673359883466679 -0.0082651356663772882621 0.1895522695283759218 0.018582516523617864657;0.026233404722042295842 -2.7699501009357541115 2.7791082959698392152 3.0674531915148350869 -0.25989572645962111253 -1.2706853375520177352 0.16937824706991211343 -2.0255767965442235123 1.1937301468536216653 0.11344705728693649238 -1.9720014508092595928 1.6121524391849040292 1.0221700861230427204 0.45918200818843524225 -1.1364020086703168211 2.0873159748854459572 -1.7445747171947685317 1.1522147016941102216 -0.2188151824986423466 0.14890252118404900816];
  
  % Layer 2
  b2 = -1.7625907466236760968;
  LW2_1 = [-0.38311629293668053586 -1.2420973593445201466 3.6709051614366363658 -6.1656516502303810867 -5.6370412218915211255 -4.7794667577872180786 4.709461721101078524 4.4311126053119158641 -1.5124550458296597277 -4.4223906348904273855];
  
  % ===== SIMULATION ========
  
  % Format Input Arguments
  isCellX = iscell(X);
  if ~isCellX, X = {X}; end;
  
  % Dimensions
  TS = size(X,2); % timesteps
  if ~isempty(X)
    Q = size(X{1},2); % samples/series
  else
    Q = 0;
  end
  
  % Allocate Outputs
  Y = cell(1,TS);
  
  % Time loop
  for ts=1:TS
  
    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1_gain,x1_step1_xoffset,x1_step1_ymin);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = logsig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Output 1
    Y{1,ts} = a2;
  end
  
  % Final Delay States
  Xf = cell(1,0);
  Af = cell(2,0);
  
  % Format Output Arguments
  if ~isCellX, Y = cell2mat(Y); end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings_gain,settings_xoffset,settings_ymin)
  y = bsxfun(@minus,x,settings_xoffset);
  y = bsxfun(@times,y,settings_gain);
  y = bsxfun(@plus,y,settings_ymin);
end

% Sigmoid Positive Transfer Function
function a = logsig_apply(n)
  a = 1 ./ (1 + exp(-n));
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end