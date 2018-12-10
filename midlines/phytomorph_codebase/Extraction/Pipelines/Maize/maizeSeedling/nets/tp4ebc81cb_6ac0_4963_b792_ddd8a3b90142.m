function [Y,Xf,Af] = tp4ebc81cb_6ac0_4963_b792_ddd8a3b90142(X,~,~)
%TP4EBC81CB_6AC0_4963_B792_DDD8A3B90142 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 05-Sep-2017 03:20:47.
% 
% [Y] = tp4ebc81cb_6ac0_4963_b792_ddd8a3b90142(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-6506.15769134605;-1969.79757456081;-1377.69545894784;-758.119043371399;-1005.44642162532;-629.490616020366;-447.378118638315;-420.527401865319;-682.765477573465;-606.219316949404;-622.812481590516;-432.343869290215;-387.57811740232;-429.96156000111;-321.30532996289;-359.072122656298;-277.037336897484;-301.07831547396;-400.794724948404;-285.384637074664];
  x1_step1_gain = [0.000244911356138662;0.000507837665104201;0.000859080031338269;0.00102992699402406;0.00116476174985821;0.00153608271588201;0.00186490882819871;0.0022956472523454;0.00140550898923531;0.00187872252206335;0.00209870879134677;0.00237421637966939;0.00249683918257885;0.00271795669482116;0.0028776513619479;0.00323536160952192;0.00314965135524182;0.00399550556362186;0.00295376821702633;0.00375051214875968];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-2.6598873625125878384;5.0055339461581311866;-0.4123293688932360701;-1.9728230352481859189;0.028729720873880384979];
  IW1_1 = [0.32705539298882552179 6.9961001926357475611 -1.0101667073150524256 4.3139581017804236041 -1.0461288834040751716 3.5841753506682576536 -1.5729018682394164674 0.4162283902040715966 -6.220391740504064515 -2.8548877712365969828 3.8267688726947879374 -3.0187574766649407643 -3.1114818234240950723 1.7041137959948671199 -1.671593273522985923 -3.1327506988148652312 0.97747882491635662561 0.79571985227910990712 4.607159929474099691 -0.97437628676543031059;-0.32259625289121407477 3.8062451952933935573 -4.1393779653439084498 -0.010217558052333082161 -0.0318666651627063785 1.029917192186549979 -1.1322392367760849563 -0.28814608940948327476 8.0968814007612710526 -1.2070859830625162434 -3.6122556019963676377 0.95283308238977271198 -1.952441531064110336 1.4766088883403702958 2.2370283846479304124 -3.1058025345534896466 1.7309245313331407257 -3.3829911220483368517 1.8672957027371028804 -0.83049249296335669879;0.46595000310905604346 5.049454981945044274 2.5069798051358858082 3.1247439481607228728 0.27924387729646887868 3.2947483970384330298 -1.1302571743735601029 0.66456029009716455214 -3.5159661956747036982 -1.3444705288956699807 0.30866957988453530559 -0.041886223383629427042 -3.3919460029090617859 4.814595885025270583 -1.0906559304842078362 -1.3212880965978364145 -1.3731148141947639996 1.6653649065608540081 1.5953103208425678528 -1.0281117087618494299;0.28360873285337473426 5.4259726769851432948 2.5698389129646490581 3.0815274951630553524 1.0598595986780701139 0.27931896003866896194 0.49773283603972978684 1.9338416648652507757 -0.78137268140861704246 0.38044657510433138592 -1.7920364833373723101 0.83772172800421051964 -0.16052941451737468004 0.39502986746727269951 0.80472306451794062276 -0.68625235638506321134 0.97491630060329592666 -0.92707731371672064125 -1.0389664745095368303 1.2106803679192335021;1.5263520015506542027 -0.0076989450752066703257 2.2097612520750695708 -1.1899771767339644146 -0.90038462512016148587 1.4818841293895110844 -0.74077145768298657824 5.0793175230449385538 -4.6430154944700703368 2.9799388689444548106 -1.6623662748587395566 2.1207491615558620168 0.78956450403072153676 5.4782175858066741725 -1.4305094192576464085 1.5538557414836822179 4.5419229824862883405 -0.51844694944405322889 -3.3358830303732132094 1.3710674489412237875];
  
  % Layer 2
  b2 = -0.47755319725191885949;
  LW2_1 = [12.486848800759583256 10.060573127166373197 7.1379414733796640036 5.2615060383389060661 8.6663260598171412141];
  
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
