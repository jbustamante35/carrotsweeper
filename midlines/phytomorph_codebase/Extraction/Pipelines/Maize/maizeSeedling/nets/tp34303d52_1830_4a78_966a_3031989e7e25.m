function [Y,Xf,Af] = tp34303d52_1830_4a78_966a_3031989e7e25(X,~,~)
%TP34303D52_1830_4A78_966A_3031989E7E25 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 11:50:21.
% 
% [Y] = tp34303d52_1830_4a78_966a_3031989e7e25(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020144;-780.471896742374;-460.205555216351;-234.31060487846;-540.291178237134;-572.142715533789;-192.819908210793;-263.800490630851;-210.82077992178;-212.78401520665;-343.803894881911;-169.741999038736;-193.663740648594;-194.470619647334;-225.820325675025;-213.71944994576;-260.914631799831;-135.174215342513];
  x1_step1_gain = [0.00034254044244578;0.00027748074658588;0.000806013241167018;0.00156677380853905;0.00168393926281428;0.00269046527889069;0.0019808499204541;0.00236361816544284;0.00245784730269614;0.00195235680232595;0.00282106579415897;0.0042465807948296;0.0031726176926076;0.00555458687461902;0.00429724736478084;0.00484270636687469;0.00333931608099043;0.00598913897545094;0.00384783002604041;0.00642488839732218];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-0.064996128505262151798;-0.42303479124705795078;-4.8722645778302835495;2.5306409770577560359;0.22527331983682272187];
  IW1_1 = [0.35903075965715619544 2.8461927340419359567 -2.3780458442150336218 -2.1568924757820941984 -0.23264373937574328588 1.1544154380500954993 -0.42296232247244169855 1.5427123369150896881 -0.08766280253342462192 -1.7698764188094591177 1.4893793902434606657 1.203039654030077843 0.024409948652138153263 0.10671829871233357023 -1.283359631369712428 1.5722755214568664428 0.52802771269248438557 -0.58430778764113810997 -1.8276585257614286206 -0.1579226792546359337;0.082453291517406698818 1.6809735339911549001 0.1092742990759864169 1.1832002102797092302 2.187998812730116871 1.1778078774872082768 -0.066489577145599970676 -0.89976388506967330816 -2.5253261408317331238 -1.7662436250130040261 1.5703995363325589096 -0.57073167597158380726 -0.81308300805882693929 -1.2245287497422387979 0.78334988219690682154 -0.48646584817953203439 -1.3081836152922678984 0.36762761362452850422 -1.7871358741633751599 -0.45980822176131902479;-1.4550623795092689772 -2.6022215663219459003 -1.2680156700034224126 -0.78552860391104473159 -2.3665761982753474868 -1.0060773581276067912 -0.91548489028888568875 2.0090516556468571707 -3.0499572655402724664 -5.5015536743858746149 -0.84634437749268187634 1.3316756590807401484 -1.6620293786045070572 0.49680387541532106122 1.6481781433071731602 2.0196967383772728155 1.9373551674496767294 -0.2493633037891541393 -4.9751278442602293595 -1.4208278063261359225;3.883284348280203524 4.6424627294761116048 -0.37719957987467889193 0.42308558416423591719 2.3864352772176395945 0.80870251960065819841 -0.16310254366108462021 -2.472556042990819325 1.3967933529214673349 -2.242170756736043824 2.9347309397125278529 1.6676597830365631481 1.671642115976417875 0.36017678131628028115 -0.49042488837008535452 -0.025777182435356975343 1.9424540304575605987 0.33585480888857294213 -1.2998127274290571798 0.76708970183566949075;1.8155511359150160455 1.2476972619958341948 0.35206630145239581342 -2.1316319045769764884 -1.2867032306594481383 -0.98414180100560222364 -1.5390787362511000325 0.25554070278879403855 -1.1229417632844469477 1.4965126717919463584 0.77682760366333958313 1.045234011087422088 -0.95029136280461212749 0.99862734958409282715 -2.2822535326771675734 0.65540272205827265761 -3.5601821963223305723 -1.7644333836372534385 0.49673924540676589467 -1.530991880999998056];
  
  % Layer 2
  b2 = -5.6182023536077165815;
  LW2_1 = [3.3241648315513838696 5.2861590172160948242 -5.5073857175467315628 4.0994934669872309385 3.5824071899003153696];
  
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
