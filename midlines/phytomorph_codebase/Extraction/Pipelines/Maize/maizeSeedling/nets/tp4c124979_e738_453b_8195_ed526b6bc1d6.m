function [Y,Xf,Af] = tp4c124979_e738_453b_8195_ed526b6bc1d6(X,~,~)
%TP4C124979_E738_453B_8195_ED526B6BC1D6 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 12:42:53.
% 
% [Y] = tp4c124979_e738_453b_8195_ed526b6bc1d6(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020143;-780.471896742375;-460.20555521635;-234.310604878461;-540.291178237133;-572.142715533788;-192.819908210793;-263.800490630849;-210.820779921783;-258.18312203353;-343.803894881911;-169.741999038738;-193.663740648593;-194.470619647337;-225.820325675015;-213.719449945763;-258.8588473988;-176.115238468632];
  x1_step1_gain = [0.000342540442445779;0.00027748074658588;0.000806013241167018;0.00156677380853904;0.00168393926281428;0.00269046527889069;0.00198084992045411;0.00236361816544284;0.00245784730269614;0.00195235680232595;0.00282106579415897;0.00424658079482955;0.00317261769260759;0.00555458687461902;0.00429724736478088;0.00484270636687457;0.00333931608099054;0.00598913897545092;0.00384783002604046;0.00642488839732234];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-3.1922806891720840738;2.2106925424719596762;0.93233857882639326053;1.4182485012292669335;-5.6500692492849697857];
  IW1_1 = [4.1807136922669059231 2.3692299766072997258 1.3191044318019791159 -1.2360396173618100057 -0.94236204504956622152 0.50263962336241274631 -0.95411581358803010389 -0.67420974946738887468 -1.090486398116837119 -3.9730709917637687667 1.4848713060825911025 -1.9249912093760013221 -0.10306227665672162974 1.044876690300760691 -0.1695752511022352027 0.30777382831706528998 -1.3924495337726803257 0.39360615773782198223 2.7590197242655771781 -0.66571173937197380965;-1.1877159328202786526 5.2941739040940811734 -0.44488824528417125759 1.0076731375295859827 3.1172019766887921755 1.3941916890339858259 -0.92069312570428840647 -3.3285119092329895629 0.46514237621781040488 -0.86307809397266366247 3.2244689408704991784 -0.6693332488648335854 1.1479908526598872065 -1.528494706065922637 -0.15105679266146787576 -2.4175559769289134415 1.6467179802926907861 1.2734501556718680337 1.0378340615081669274 -0.88253251840228552627;-2.6342478456181903645 -3.2124594328602080751 1.5947765969395337748 3.4488858363575212174 0.58518430555334077781 0.14964212601194315511 0.8866120176808529596 -1.0161792930062489848 1.3050140458926722253 -2.126849291207503434 0.74246123127108398787 0.51326703635031478701 0.68326977350745565776 -0.1598544887786996771 2.7847903095392130268 -0.67143195705124836259 3.2851353109273837561 2.0290867738560125488 0.96029120922434230145 0.79078245501939437023;-2.1256349058781336403 -3.1047175049674979874 0.92575022511293980632 0.47659809044536027844 -1.6584290961493757255 -0.10877121601851508226 0.73569316747624291519 -0.33054171792303133071 1.8729604332265306077 0.34251147794147912551 -0.66306548460275516188 -0.28780874283801644031 1.3455882887110086621 -0.51207402924099199026 -1.7414304967494085297 0.75588038235300281098 3.7347461111262809652 0.67537251085750582558 0.25670122028883263221 -0.1589987823544766754;-3.5825231214911466004 -4.1514805425258591853 -0.73549607765725277808 -1.5295944740887743585 -2.5693870782484284909 -1.1083662620120464659 -0.18532585899446635258 3.7811340920339451266 -4.301582532858905239 -2.0559724507877450428 -1.333304448972249423 0.37001382606820470977 -2.338133424895568524 -0.1722467084739592591 2.1983213935378471682 0.51061936066311908178 -0.21819562470233766494 -0.53427559215036390938 2.7364920442558569036 0.7313163474761542604];
  
  % Layer 2
  b2 = -2.7123205765876234885;
  LW2_1 = [6.6020608465047283531 8.5711783581084954875 -6.5110029238294044873 -5.1566766519960784976 -9.6467754369288662986];
  
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