function [Y,Xf,Af] = tp057ffe7b_f706_4794_9500_9308a88378de(X,~,~)
%TP057FFE7B_F706_4794_9500_9308A88378DE neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 04:36:20.
% 
% [Y] = tp057ffe7b_f706_4794_9500_9308a88378de(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-3755.36764584613;-549.061967680733;-873.097427232646;-1646.05375716082;-423.333970720151;-622.824124154291;-251.161393954612;-450.689400777678;-342.308065266756;-280.547165434121;-739.082929249962;-243.822193480699;-418.782145732691;-227.785437026336;-141.910346711069;-156.623146275712;-201.34735048819;-147.229750402556;-715.253580951428;-185.847378023652];
  x1_step1_gain = [0.00043525835929838;0.000254266760874673;0.000347353003038839;0.000777513986910233;0.00172505617054034;0.00206247464333953;0.00302249934878052;0.0032309825223835;0.00220614635584285;0.00242812082773061;0.00231897377235304;0.0040749586772638;0.00341201206501295;0.00551378415971665;0.00648970220781285;0.00620026189110804;0.0063360245550118;0.00774864116842834;0.00243421857993076;0.00313801600663603];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.7761365166419085249;0.76331778554252460722;-0.15703190832341926142;0.67110462624354894068;-1.0113791837632968296];
  IW1_1 = [-0.020251036452187155812 0.26002448015311346863 -0.14920534284369152411 -0.28818548206715510762 0.58839856156180392333 -0.57953105270930205251 0.42134543158047566003 -0.30255044933314956257 -0.20432907557610693883 -0.094110431174067199489 0.14328198530770025476 0.18093203759476034209 0.17213758329140832481 -0.24062861965371845629 -0.18482000766810205428 0.54294585286590246742 0.38548757046831572115 -0.17241702374267189879 0.1281278127987900417 0.55351058904274708095;-1.7944895152074460398 -0.047123721240247738729 0.6717553842873899006 0.5259475966636394384 -0.096666860851147409783 0.56564007974344132634 0.20315337158223659153 -0.15948593503106761604 -0.17010163006143719233 -0.49250115678344902337 0.28459879546052585519 0.36742926748497406475 0.087590511207620944112 -0.25159643144124127501 -0.45205200182842542045 -0.29951537955692153492 -0.047910243195857710186 -0.30724806557235473381 0.50985306543489750108 -0.11547829283011702695;-2.3969422605125432213 -1.5520817835798157436 -0.45655512939367054459 1.5971229745432926617 0.28711276457277679075 0.051606312648067961812 -0.61217634180777791286 0.89438865177449466959 0.0036392212316548161986 0.37821861176295945572 -0.023304707281682872366 -0.44192132118301707955 0.10505641293306300366 -0.14128946706378681109 0.20062936930493374676 0.18921795115268830179 -0.49986151242026755082 0.30566150792524737767 0.22251117749820770175 0.45907916816368488089;-1.3946192136567714392 0.39635487064332114437 -0.051341963503906115629 0.27727050022792698591 -0.14955622921123684632 -0.15223931176305752211 0.069606518382475002227 -0.095629070163576995744 0.10753829181904184242 -0.30688291369980180701 0.38954329906222096014 0.48042368202924801146 0.05914200549322814332 0.39167735532816794208 0.13654541401131678469 0.036640718348260277459 -0.50239299049865626667 0.40994010295261934518 0.44235105263575097423 0.48421671335395610747;-3.0167439159252276148 -1.8114170607158499582 -0.62663548734248164518 2.1497268906312649861 -0.28193732791563025364 0.73052450761342924146 -0.11029436405409633049 0.63494631839974313703 -0.28534726169891327929 0.12782701376982463781 0.5551797244105844964 -0.39195840130999859419 0.31934114340523211961 0.14125399349207307109 -0.14524861116699061547 -1.0272987693793924624 0.45026861932108019904 -0.26158403744413694092 0.26959563572448652646 0.17938139271371328354];
  
  % Layer 2
  b2 = -0.31488639972156690083;
  LW2_1 = [1.3178303291118482754 1.7201690581247524037 4.0342386661643896417 1.6244984670181661368 5.1812737131186672457];
  
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
