function [Y,Xf,Af] = tp102246a4_f937_4747_a7b8_846f6ecfe8c4(X,~,~)
%TP102246A4_F937_4747_A7B8_846F6ECFE8C4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 13:48:20.
% 
% [Y] = tp102246a4_f937_4747_a7b8_846f6ecfe8c4(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923585;-1028.71862143182;-451.710971542892;-692.538815519986;-441.04379323221;-783.667639824722;-450.02627128993;-656.555828069294;-380.160028036838;-183.58947139104;-268.743398235208;-384.146123020656;-202.979007080836;-313.061223281942;-302.321350914787;-296.41990859127;-275.821003918814;-174.982344511311;-176.062817024264];
  x1_step1_gain = [0.000192953812201724;0.000635170113819406;0.0010737117105794;0.000992967389065702;0.00140247416929524;0.00179426803546527;0.00152094087103249;0.00262475201653713;0.00199415494766508;0.00251075910029926;0.00454433874041192;0.00378662000738747;0.00302237737948978;0.00448879604138558;0.00366080733460348;0.00363057230533701;0.00381924645546724;0.00427537582274535;0.005921205517812;0.00602422138644422];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.4770671508018757834;-1.7296347882079341485;2.2858052596027760295;-2.8571641798541280721;1.3370037709628304068];
  IW1_1 = [0.32883443252161426162 -2.3012208941193539324 -0.4754297735397867708 0.13872886270941131581 -0.91754450734620374153 0.22030565724708212061 -0.081200311789518420102 0.38262280969132383968 -0.28303720828258954745 -0.13509760444038007576 0.50678130144589084161 -0.71874136890729123373 0.23874429308161107532 -0.54497962048074688468 0.41566127994053320149 0.098591363897441275377 -0.43043575970819436538 0.29415284815673364838 0.2878270903225138011 -0.11623135303848323996;0.82229293191425167908 -4.9352296376266151512 -2.8166422097396801583 -3.5113226420051439902 2.3913448209472472072 0.87776803719763896172 -1.4203070522110836649 1.1452219570006760829 -3.3366791491754344001 1.36701370141787093 0.56952385469921884553 0.74428638492821397321 -0.054284101133596079913 -0.15627983915393112158 1.6981567365828873495 -1.26581824616401728 -0.7169393588381768323 -1.5100893600140750106 -0.20107287866811457411 0.15271936139571490321;1.3938809448208306918 -1.5281629755971817275 -3.7582013095978528128 -0.24723745321920939655 5.6999794756158221887 0.86433003602863023307 -3.0147747740193278077 1.2099124159171266157 -4.900752070312118569 2.4821655164086289247 5.0685524162273569715 3.4644357836167225528 -4.4471969596705509531 2.6302565840514700035 -1.0816772908616849591 -0.8535543399405763898 0.80279031088556407347 -0.73522725193758331486 0.15122781556236578049 -1.0458507338421529909;0.14724012499702929246 -1.9289139887717507982 0.061411710102226096941 -0.42062576639459064998 0.17904545582701464346 -0.7730090058720511248 -1.0730957550220840968 0.031995247021593745118 1.0284526682662973318 -0.65247428173501076376 -1.0110123319115269513 -0.43023852791451555611 0.6482323629665152076 -0.39039540005186962324 0.68962325655326350837 0.37719640444288798919 -0.092294148935446962412 -0.96974110878348129816 0.59217888832923915032 0.15832561905583936168;-1.3654810202643363759 1.4478237040680990155 -1.8573635439174525086 -1.7552358688517932794 1.7110418827120543739 0.19236264745707523649 1.3761575203691640557 -1.5015008869997699126 -2.2254787820791741559 1.4716478929370586659 1.8664664063009486217 1.0629544885507342666 -1.5757862637301744968 -0.29980044679503486948 0.48351083593512050873 0.92751057318365182081 0.55839592371293478301 1.3029231194186297049 -0.11819986200601524662 -1.1612238099487566245];
  
  % Layer 2
  b2 = -4.4626328740479115353;
  LW2_1 = [-1.7859113562604329584 9.3152172499704004593 -9.6281583062480393664 2.9559159396436935907 -2.5133981066789186087];
  
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