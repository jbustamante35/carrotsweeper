function [Y,Xf,Af] = tpe0b56f95_a10b_4efc_9748_1c972c12e99a(X,~,~)
%TPE0B56F95_A10B_4EFC_9748_1C972C12E99A neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 04:19:27.
% 
% [Y] = tpe0b56f95_a10b_4efc_9748_1c972c12e99a(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-6506.15769134605;-1968.46879957595;-1377.69545894784;-758.119043371399;-1005.44642162532;-629.490616020366;-447.378118638317;-450.68658711172;-682.765477573453;-458.333844314355;-622.812481590516;-410.039334332747;-387.578117402332;-429.961560001115;-321.305329962886;-359.072122656305;-277.037336897484;-301.078315473927;-276.306473605863;-247.875867270408];
  x1_step1_gain = [0.000244911356138662;0.000507837665104201;0.00085908003133827;0.00102992699402406;0.00116476174985821;0.001536082715882;0.00186490882819871;0.00229564725234536;0.00140550898923534;0.00187872252206336;0.00209870879134679;0.00237421637966938;0.00249683918257879;0.00271795669482113;0.00287765136194787;0.00323536160952178;0.00314965135524186;0.00399550556362225;0.00295376821702617;0.00375051214875961];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [0.57954289003952974557;-0.25388905969801817131;0.65762748875697163431;4.3341239612310280549;1.8722915713679415539];
  IW1_1 = [-0.89162470191103326833 1.7648535979203752166 -1.1354238189176166784 0.62123072523175226767 0.55617961352970191058 -1.7500312475219279573 1.3185995211977805663 2.7837460940035123613 3.7625622439937060904 -0.17559868876443285313 2.5363717198919717433 0.62279611502649256938 1.4068057584454347619 -3.9974597161989637861 2.2633914520095559375 -0.70911205149461431496 -0.7924257045370933783 -0.95383146149136266434 1.7309070829721138196 -2.1133364241747836232;-1.5414951142877253432 -0.32180907031324912371 -0.94236130208443702472 -0.97412599252403708761 1.4685180108810511612 -0.70084793176282778049 -0.0081659263462928489974 1.3834029807836005777 1.6298008048512295787 1.000141601345491793 -1.1365389794730449413 -1.125645041071248853 0.7090215684288938558 -2.5324776765553860258 -0.61028240038855885441 -1.2923074583222198264 1.2755053641998395086 -0.4395170631647333237 0.2174150211134296129 0.12980047399099772787;-1.0071996163700664972 -5.9639102756977981912 -1.4637892942035561905 3.5683499595625645107 3.0304591213656006943 1.3822805401211026854 -2.3390026591102541786 -2.6250701549676080937 2.716199585120137705 2.2636178543428422394 1.5720624431707050928 -0.70433987232426953362 1.297799978781651653 -2.0343728800344615948 0.98255677330915258683 -4.6723054937203585268 3.6048181854323417816 -0.98854149152459325212 -5.3496005795110415093 -0.75206855190280497148;0.48471864072064013707 -5.3262290591282894781 -4.3244791548634013978 2.7462125460080173411 -0.09906857904911514634 -1.2002683794773483239 2.7763289389939167506 1.5218661932103054735 7.2508145441187492608 2.4437000030086326596 -1.9480887525203425259 0.023555692246005611318 -2.9115331037736589259 3.7186222561922583907 -0.3846002176488986235 -3.2772042366642533118 2.1920914006407601349 -3.0500450044274254324 0.44631554812190615822 -0.82303409337340593943;-1.7150830159522105767 3.1137030540270105305 1.7153778471209932821 -1.4487881981827825939 1.1605513169165502241 -1.4783126199175078419 2.1083262836966949472 -0.39280824872816155757 3.8497517655627100375 0.26559001766486783902 -0.19033275944744090791 -1.645086063604343396 0.87244221190590687254 -2.5948664771714446964 -0.18138999017301288275 -0.30200207219490865507 -0.65917386335166505695 1.2124084828287522431 1.1006460778774154186 -2.4111062853860194544];
  
  % Layer 2
  b2 = -2.2010090773277868337;
  LW2_1 = [-9.465595675170494161 -5.5131745446825872747 8.2209405885307038631 8.7941448155615216109 -8.4332466640725467499];
  
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