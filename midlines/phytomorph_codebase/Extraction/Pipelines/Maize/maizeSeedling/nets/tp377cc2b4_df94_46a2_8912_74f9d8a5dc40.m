function [Y,Xf,Af] = tp377cc2b4_df94_46a2_8912_74f9d8a5dc40(X,~,~)
%TP377CC2B4_DF94_46A2_8912_74F9D8A5DC40 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 18:12:53.
% 
% [Y] = tp377cc2b4_df94_46a2_8912_74f9d8a5dc40(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5342.40771051541;-3062.40716175129;-1831.68761111049;-751.778745457493;-962.879682501258;-561.234234396122;-567.818719457228;-566.708636588772;-460.813029179954;-650.116847383037;-365.703888775152;-497.927214156977;-445.050021503471;-362.02756517138;-290.108959264264;-228.419697429558;-268.266206223882;-201.923452595307;-223.562425332956;-192.300211130342];
  x1_step1_gain = [0.000187746732766858;0.000517587446146427;0.00076092913976903;0.00107683880161886;0.00113410115771159;0.00154637953830619;0.00191937835803698;0.00192268942981793;0.00211835684506625;0.00196958059572072;0.0026895539866255;0.00257377430529597;0.00284234016821654;0.00279908218986806;0.00409640929201521;0.00345584818973336;0.00422505290367564;0.00549698021899224;0.00471432618261879;0.0047687856727025];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.8862935923658235016;0.92233611487233924109;-5.5190259028811450293;-5.0276402237052613131;1.3511349229399116645];
  IW1_1 = [0.091840436069403300068 -1.3453605359750000492 -0.43577011104435181554 0.31565545192903082805 -0.2651708699852020068 0.14232653399781955428 -0.094091217604830551524 -1.2429876875741387998 0.41602331966079392211 -0.065547576796555165335 -0.32719355445815689887 -0.39503088040014366911 0.44038766429133091584 -1.0012036407433428931 -1.0951453638273687119 -0.76357166159908951286 0.20554793037741458939 -0.16498740644494977459 0.39999610607515184713 0.35741679374292206761;2.7437293546043672876 -4.1437708933467778394 -0.035676857461289498974 -1.5233644593637150066 -2.4979103555330506303 -0.041114569038715999272 2.4979371535715202413 -1.8809782155365386558 -0.12514430579419624778 0.449284985301904638 -0.1225569012891088716 0.1157860523199684688 0.17368383434791712094 -0.51245629986422269742 -0.63567991318669536671 -0.71068009707703549971 -0.25008610318735685629 1.1281141511013985657 1.3841290131338863656 -0.085539284040282048682;0.28812191871528053388 6.2328322106480529285 3.0125699333049080764 1.4207543274124274646 1.1503894601622932647 -2.3131263214760329205 -0.55007583389422776143 1.8333550908375135791 1.2114847774677290904 -0.90267338142714248495 -2.4322191278647111723 -2.1241534146711029685 -0.57460897110019937273 3.1447788147401651848 -0.2193419854063466401 -0.73451955710129079868 -1.4099329072918993511 -0.14124760962257337793 -1.8638847052198914689 2.7508853203353176298;5.1218428266604174581 7.0802386117775029106 2.0019314353169423981 -0.49132627010611917839 -0.62629334864070007605 -4.699501731006190397 2.9258224579787777131 3.0105040547056374223 -0.022198520389327021374 -1.8974910510667883923 -2.5783269138262054909 -0.48681523270021370697 -0.41161507444488271146 4.3876286809519493559 0.76669378616845729191 -0.35359807651810593043 -1.6715636035925218739 -0.14837959817331664025 -2.8861576092326242104 2.5908940795556381076;0.31222175005606006559 0.54653205166721041675 0.21301107524103099489 0.48059024902174990457 -0.39873845894498544107 -0.58591138509415197078 -0.13541641874796006562 -0.39581149345506955184 -0.32980290125676536839 -0.40236249843839783091 0.43754345699302954475 0.16686736341645230786 0.37784828744962895009 -0.52711687646645732919 -0.46720436875703991175 0.004617331238185409456 -0.028889337727143523521 0.20772901601655224813 -0.35320915214220577338 -0.21338905979219410791];
  
  % Layer 2
  b2 = -2.3390514829709267985;
  LW2_1 = [-2.9748559835193542078 -5.6136451540165897001 10.134251758430314894 11.704385816213916272 -1.1595881453497780544];
  
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
