function [Y,Xf,Af] = tp70432afc_48c5_4253_a82a_8c61923c422f(X,~,~)
%TP70432AFC_48C5_4253_A82A_8C61923C422F neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 16:42:31.
% 
% [Y] = tp70432afc_48c5_4253_a82a_8c61923c422f(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923584;-1028.71862143181;-451.710971542891;-692.538815519993;-441.043793232207;-783.667639824726;-311.950474474655;-656.555828069289;-416.411813428014;-256.51856505755;-268.74339823523;-384.14612302066;-202.979007080799;-233.266353450959;-248.555985162485;-227.243600250016;-191.974128700863;-174.98234451128;-176.062817024317];
  x1_step1_gain = [0.000192953812201724;0.000635170113819408;0.0010737117105794;0.000992967389065709;0.00140247416929523;0.00179426803546528;0.00152094087103249;0.00262475201653714;0.00199415494766507;0.00251075910029921;0.00454433874041207;0.00378662000738733;0.00302237737948974;0.00448879604138593;0.00366080733460337;0.00363057230533686;0.00381924645546832;0.00427537582274492;0.0059212055178125;0.00602422138644343];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-4.8077218790205495935;-0.81697960142450332111;-4.122066878063090023;2.5405907619693861577;-2.2470163884355311978];
  IW1_1 = [-2.0482889931890704283 -4.1130039179907429414 1.8415219479856150819 -2.9183845742931744383 -1.8166990199081893831 -1.1514059421285620655 -0.88460050140175527922 -1.1678431027491038563 1.5353258974596881714 2.1719938441658119821 0.53302959662629711968 -0.085441846953706448997 1.2637475085706599831 -0.017916147096144058315 0.66169099981196988214 -0.63875810535063282103 -0.16494454456185042579 2.2439642841847331489 0.81340867815039785693 0.88858647907019361334;6.5712573523070849646 -0.6220863768042067532 -1.1937038749456621023 -1.2905993471011236728 -0.28690336853435294273 2.0764784937455731395 -0.92435296191814897426 0.10929961771517639202 -1.4335194992707469552 0.20372637995316086879 -3.5619701652629420785 1.6259703883502227839 -4.5004222105775149032 1.7124208035707575437 -0.30544018663774002986 0.053335451783630347111 0.35829376093026071848 -0.30782657924345718792 0.58252694434854646044 -0.7059136243415350398;0.66233525923490832454 -2.1401824713371491171 0.34603385470314884653 -0.36296264862483468283 -0.29475225773295182785 0.25111641304449772161 -1.902508104991801563 0.40412962730875401363 0.34374638838299864663 1.421480864640437547 0.37873744565456607392 -0.13872512284949817252 1.5309535582202444814 -0.41305187046794267225 0.11293644467421477817 -0.20136789722244333167 -0.14620468372131906198 1.4792938007885474416 -0.20492070541187970556 0.75488415394202379449;5.4683313936901800645 -1.2726860367826045017 -2.9367513778190450502 -0.53992727980804477994 0.93315733854264426572 -0.23894323155819136018 0.18275271996361086568 -0.41701115608247957267 -1.6587671074805474536 -1.518344746041512705 -0.94553000815757548558 0.83964014664542296895 -1.1684200560911113254 -0.38595886349883445909 -1.130243633529409264 0.74744222842648400107 1.0613312188079266907 1.3214831211802824384 0.49327767876408207259 0.33654694776854804372;-3.9164548190857879639 -4.0556433908554865297 -0.11448179228941002261 -2.2695947006803529611 2.0981652555482210865 0.41036289646556611821 -0.64893442782427290805 -1.3426817715225607408 -2.5764705426950929024 -0.77032298555459188805 -0.87460618567834114145 1.4650435104653558582 -1.7798674914384424728 1.6576898397488424308 -1.3558435127305588086 0.18271820906768462822 0.55853341357772734455 1.0278677313295923046 0.593365205758598524 -1.627756136973801171];
  
  % Layer 2
  b2 = -3.8579088597350890133;
  LW2_1 = [5.1312295796261571112 -3.2531328203680520694 5.2321386097100912593 6.6770553631752269652 -4.2394792282553313001];
  
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