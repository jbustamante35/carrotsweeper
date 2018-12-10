function [Y,Xf,Af] = tp1e5525e1_a801_461d_9247_3a540fb00d4b(X,~,~)
%TP1E5525E1_A801_461D_9247_3A540FB00D4B neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 01-Sep-2017 11:43:11.
% 
% [Y] = tp1e5525e1_a801_461d_9247_3a540fb00d4b(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923584;-1028.71862143182;-451.710971542891;-692.538815519991;-441.043793232208;-783.667639824725;-450.02627128995;-656.555828069294;-380.160028036838;-183.589471391039;-268.743398235219;-384.146123020663;-202.979007080814;-233.266353450973;-248.555985162466;-227.24360025002;-191.974128700923;-174.982344511177;-176.062817024381];
  x1_step1_gain = [0.000192953812201724;0.000635170113819407;0.00107371171057939;0.000992967389065705;0.00140247416929523;0.00179426803546528;0.00152094087103248;0.00262475201653707;0.00199415494766512;0.00251075910029925;0.00454433874041203;0.00378662000738738;0.0030223773794898;0.00448879604138527;0.00366080733460352;0.00363057230533714;0.00381924645546787;0.00427537582274354;0.00592120551781401;0.00602422138644287];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-2.9529198724058045578;4.3228622856450424905;-0.10552286588436896042;-0.53867901205649293406;-2.6445773953381448429];
  IW1_1 = [-3.6965200980762169536 -1.1913512166359296884 2.6531187487968455052 -0.24597672605974607585 -4.1678321015245956005 0.18529866879117606882 2.8948360914592905502 -1.3055798789917476554 2.8382110141125278524 -0.71058910216831205897 -4.1426122972443382153 -4.4341122696246344503 4.5493214147014411353 -1.605192188182358759 1.2395510414721977455 0.17278335427924865209 0.69340516505282223658 -1.2747054500255767806 0.33547259303020998678 1.3360557799275427637;-1.3645536942376572664 2.5188613478382300137 0.85572948862106823498 2.0058206126592734719 0.6449206388521319111 1.0637058127310501376 1.4048493400404551412 -1.4161480132164354639 -1.4842116027745746454 2.4536581967280928218 2.106476344164792458 -0.85952431985862443042 -2.364577195254735642 0.4671720701107763829 1.0916647980917011385 0.90581644489966850475 0.44952436695897007857 -1.7600143681066620083 -0.22894759559606298049 -1.5645126970186453708;2.6374542473495763417 -2.5181176505953928313 -2.3644785552288425379 0.83315972384664027661 1.5206366545255061151 1.168492043893475385 -1.164462706415306803 1.2275407741116279858 -3.1312881804582310963 1.1964064868864672242 1.4078871994401127576 1.2112062055057148946 -0.99211190127481330148 1.0388574458771784315 0.78671533872645571517 -0.17086897092180544644 0.18067716277844952932 0.76652830697976148677 0.35080297043537372881 -1.0505250205675149555;2.1990984381749547616 -3.1890415095145714375 -1.8115052708973293782 -2.850063939283355996 2.9276461500275634542 0.40354219083484349495 -1.5662416725860910383 2.3499969763986379157 -3.2161958295949601982 1.0818634082454292589 1.0286364114698707262 0.85874653285840152339 0.27335094980932717057 -0.38509689422603254894 -0.33497013711239448375 0.50218826852069886613 0.77395508017743097628 1.3398052465526151256 0.34790187085653256904 -0.32059335253255366593;0.20531726094647095859 -0.28872303080588035407 -0.22282419485800597503 0.39979911901167219934 0.087512376075666448161 0.15045608566514431792 -0.13579762519930532738 -0.41711181639383482223 -0.72263992707371538238 0.37319054817576385874 -0.12750723558381163447 -0.098959750700817972202 0.16289560772574079861 0.48396454875150923813 0.23638416060294703547 0.37866405617489373725 0.66386594292835165998 -0.17518780340685269348 -0.42566885035369605994 -0.12247089141002426294];
  
  % Layer 2
  b2 = -4.5517928832004921347;
  LW2_1 = [9.4666353904861022528 -7.2285583423566297512 -4.1190984692863228389 8.7079046197027114573 5.0598679664970553915];
  
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
