function [Y,Xf,Af] = tpcb6419ff_f07c_4807_9932_2347cd17397b(X,~,~)
%TPCB6419FF_F07C_4807_9932_2347CD17397B neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 08-Sep-2017 11:08:01.
% 
% [Y] = tpcb6419ff_f07c_4807_9932_2347cd17397b(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020144;-780.471896742374;-460.20555521635;-234.310604878461;-540.291178237133;-572.142715533789;-192.819908210792;-760.602424664802;-210.82077992178;-212.78401520665;-343.803894881913;-190.320782432691;-193.663740648594;-194.47061964734;-373.104709120359;-213.71944994575;-260.914631799829;-135.174215342501];
  x1_step1_gain = [0.00034254044244578;0.00027748074658588;0.000806013241167019;0.00156677380853905;0.00168393926281428;0.00269046527889069;0.00198084992045411;0.00236361816544284;0.00245784730269613;0.00195235680232595;0.00282106579415896;0.00424658079482956;0.00317261769260759;0.005554586874619;0.00429724736478087;0.00484270636687443;0.0033393160809905;0.00598913897545122;0.00384783002604043;0.00642488839732256];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-2.3981045881345584014;-1.0867531515303299194;-4.6779633057490928039;3.0181214876801476521;-1.4371299785171274266];
  IW1_1 = [0.87503904529601239659 -1.2859876413820963137 -0.16783943614825130863 -1.4866161031769629108 -1.3054056887538270093 0.097059848694068556751 -0.46260269964942596355 2.4304095144248529792 0.42957867911656977267 -2.4543662876936109818 -2.1473460284343883764 -1.717876265917463563 -1.0013181383808369418 -1.1491016240224012357 -2.9599092295811089137 1.4532953336119232013 2.1086878484115003651 -1.237584098554108003 1.6948169711578922048 -2.3160791000723324196;2.5258915163407129612 2.3656501116777923777 -0.37858402435127580521 -2.1435969494560191961 -0.21543319408742483967 -0.099124536872010809985 -0.37023524997570489425 -1.7443344053266274241 -1.0152591719339674192 0.12753813683187764649 -1.2015286658433348421 0.23615385424902490819 0.4118469333144031963 -0.22939976535454359441 -1.1951269951297898153 0.11310771678049233446 1.4167333838248961086 0.35676397879913013256 -0.93296343734630104993 -0.20413021110221096777;-3.8343623339167525188 -2.5958597089242640976 0.23525370937241865676 -1.1097426771351004415 -1.4642945427959079829 -1.9349629810817341458 -0.44460949003688637893 1.7213101487202415729 -0.65087998251844136899 0.87083436073637077168 -1.8332418056672088458 -1.1979960546575694202 -1.1601736615387687745 0.53572908942873087579 2.2529357404183802949 -0.82549435240289259497 -0.36409321373994568694 1.5352697155017611674 -1.5861639303720871386 -1.4642134377503257525;-0.86510225759564451931 4.1302144067437351893 -0.50753637481512003937 0.16777076080742295283 2.660771245021666509 1.3584040110073705154 1.1356000458900163963 -1.2278909141102556379 2.3235491372852878911 -4.576476189587252108 1.1514474259245823795 -2.387438585766803012 0.24327473034923458028 -0.20445121882537309332 0.43839002903117935128 -1.277068795813935731 -0.55675533092369311472 -0.17805491608659582137 2.2133548466955290657 0.71904515049907857804;0.28169115819313955695 2.4960410682497902179 0.79577782878767744368 0.083839056230180758478 1.6375651565014510691 -1.2887412019363138072 -1.0877411564439707359 1.3896699978299380529 1.8520348731717815038 0.29449459935820304635 1.9206432162855562318 -0.23845229993408692537 -0.97823402821994598355 1.0737673694587308582 -0.33986977989811778533 -0.39459201398847254527 2.7605313701316189778 0.067778816345993961323 0.39888268647596003547 -1.2877659912387422469];
  
  % Layer 2
  b2 = -3.1380173565054283458;
  LW2_1 = [-4.8898758683571621475 4.8989659977712536332 -5.7535211923095479492 3.640261506285807247 5.8783192342747776493];
  
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
