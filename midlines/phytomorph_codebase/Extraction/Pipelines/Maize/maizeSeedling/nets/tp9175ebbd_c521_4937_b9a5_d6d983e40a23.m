function [Y,Xf,Af] = tp9175ebbd_c521_4937_b9a5_d6d983e40a23(X,~,~)
%TP9175EBBD_C521_4937_B9A5_D6D983E40A23 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 31-Aug-2017 13:50:34.
% 
% [Y] = tp9175ebbd_c521_4937_b9a5_d6d983e40a23(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5432.21511079255;-863.132365464489;-1030.35626722576;-451.327248716542;-690.442124301695;-440.923810847593;-781.659803211367;-311.516810675615;-659.850047510512;-416.874768568223;-250.644035779545;-272.098158795939;-389.737207411697;-200.812705665977;-311.371005878013;-243.035648213554;-279.441974279547;-276.019527491442;-178.529888662629;-167.395512111997];
  x1_step1_gain = [0.000193639844679507;0.000635136090067475;0.00107257204462677;0.000992521373387899;0.00140098678161598;0.00179249040902968;0.00152180852073824;0.00263160812648061;0.0019868610334328;0.00250748142155721;0.00458743324706901;0.00371176793861507;0.00296778207281887;0.004489150933777;0.00365879017755123;0.00377014651183433;0.00363538636865064;0.00429258503045445;0.0058738883725327;0.00612684175841709];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-0.6856387644138317583;2.202879941396193697;-1.8695694687901618369];
  IW1_1 = [-1.900393474177088482 -1.2175506420618811276 1.8457098343915210314 2.423912834580910669 -3.732575219994203497 -0.64277197020138099326 -0.44863992823863879655 -0.6163822871554320848 0.91974804088375861078 1.9485399465038042432 3.5611226456160194331 -2.5007692003325110441 4.9440678622879339343 -0.31424067787488207948 1.1422521043498519155 -0.8313265584661871932 0.52339248603124510772 -1.2419632421514366083 0.095207468911270909384 0.38921050084363517652;-3.0237243414773988803 5.7970308931077818571 0.65142578139550044725 1.0501343792264439259 -0.64793428779842687781 -0.067229654525203394466 1.6503771877775261778 0.00047341899798388972775 2.1863346498105440396 0.31534571053166343679 0.39328871702550050049 0.51544026969755807333 -1.2051961230584053997 0.19331967825334478439 -0.55720491011021744399 -0.6413925639549624691 1.0861554561007058517 0.81977957367945963885 0.10430714695139778336 -0.22320905016705630541;0.060343480227570786556 0.086303312405268450824 0.14721401755518798149 0.60308283018702646405 0.26706554355891887198 -0.40231769028510855346 -0.53091135360716390501 0.13065162957582288072 -0.052344884668740578049 -0.45321257284982890967 -0.21543527141953905391 -0.19807153055854942458 -0.46671699053435872218 -0.48557004771477751959 -0.46258251258150079321 0.42640383265307457705 0.22883616391603572859 -0.45154941413076393752 -0.03814853752923944713 -0.10019274536027934774];
  
  % Layer 2
  b2 = -3.6365332053968657533;
  LW2_1 = [9.4457429706098512412 -8.6841799295427524186 3.610195550899803596];
  
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
