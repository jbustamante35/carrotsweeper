function [Y,Xf,Af] = tpf1b0f249_fbee_4d12_83c9_8134f52996c4(X,~,~)
%TPF1B0F249_FBEE_4D12_83C9_8134F52996C4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 31-Aug-2017 08:28:22.
% 
% [Y] = tpf1b0f249_fbee_4d12_83c9_8134f52996c4(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-4896.23771001091;-863.132365464488;-1030.35626722576;-451.327248716536;-690.442124301697;-440.923810847591;-781.659803211379;-448.474762490208;-659.850047510504;-380.738319531136;-185.329610546061;-272.098158795915;-284.166721828048;-200.812705665954;-235.25777110712;-243.035648213602;-270.705877196786;-189.90018890573;-161.96006844184;-159.036943443448];
  x1_step1_gain = [0.000193639844679507;0.000635136090067476;0.00107257204462677;0.000992521373387902;0.00140098678161599;0.00179249040902968;0.00152180852073822;0.00263160812648067;0.00198686103343278;0.00250748142155717;0.00458743324706882;0.00371176793861535;0.00296778207281901;0.00448915093377834;0.00365879017755116;0.00377014651183343;0.00363538636865051;0.00429258503045364;0.00587388837253377;0.00612684175841704];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [0.5014258197819565499;-0.073878685454639395069;3.5602997249692585058;-1.6771093068209927957;1.8289074099050615096];
  IW1_1 = [0.68368591827943603345 2.46670947627567827 1.8550083144021984882 1.2980853831823653266 -4.8328235903756144864 0.32323363168044516591 3.8928119726534933243 -1.7604468772717882796 3.0876531595426079946 -1.2782376723842072685 -2.9367496255805489724 -2.4852298331061146719 -3.0409974537106077008 -1.1353745478369376976 -0.40239126128899793677 0.18835325585690620986 0.79361150483744724582 -2.0720932018008637598 0.47158696781762415862 0.16886689610988506183;0.96121714880655972912 0.81891464020225568365 -0.93720422027630601036 2.0773850860378351157 -0.2633082843971633058 -0.63314471840078501153 2.1940622436548866681 -1.4701994185025004658 -0.43005752156237375683 -0.40689225944705537108 -1.2620500483886918364 -1.3233518665460579644 -2.1610956545254871308 -1.105519125934072644 -0.39211742699054474848 0.2796233571075557589 0.19740216226488777052 -0.9073562422776350056 0.38351296903494758173 -0.66390759630371487532;0.45757808914200526207 5.8677208959447213488 0.82593890826527116555 3.0849511049294586762 -0.5777191379646006153 0.16711806363344511306 1.8313475696220880007 -1.4353977905992325592 0.71627718080755420349 -0.092510173056384675805 -0.64512827533262462332 -0.58398703844274402197 -0.070372629414207271692 1.4035476438386342846 0.70211741693021334498 -0.18430487887278002024 0.48873782847297142506 -1.8560894990366265667 0.97058154921080785904 0.49387428823841350578;4.2250417196186829116 -2.257412438137940125 0.27439065214007213322 0.46795465775463079838 -0.34054103516241662764 -0.33325222223687184231 -1.2735130422695046803 0.31079105724837635671 -0.25232898654276408879 0.49680633212255048692 -0.77289358580204292615 -0.24380859797369500153 -0.36242458208384620288 0.25591909183379829562 -1.3794881178811373434 -0.6260149450864079812 -0.044641516797856675702 -0.86181796175300418028 0.13760237486075968527 1.0522552698814633398;-0.17842422169897334538 -0.78418030023676588147 -0.52228929746692598357 -1.5094865599393261224 -1.4193711768341674695 0.61759772361970666221 -0.091195366537967240128 -0.62033819331678374809 -0.45947952024024113493 -0.71540783849247047144 1.0917331384432626962 -0.52303783389535329906 -0.66127175323110698457 -0.23172913188445648514 0.092795392449609478969 0.76855708931940691642 -0.2197677101600891969 -1.0667563383108313868 0.22155148658783407889 0.49169028835423139734];
  
  % Layer 2
  b2 = -3.9385166320790023775;
  LW2_1 = [8.2153948729162760145 3.5367485319289566625 -8.0903166463831741595 -3.976551151809331941 -4.3167282069557231239];
  
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