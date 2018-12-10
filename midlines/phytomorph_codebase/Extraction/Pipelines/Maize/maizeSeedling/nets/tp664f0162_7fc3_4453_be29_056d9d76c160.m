function [Y,Xf,Af] = tp664f0162_7fc3_4453_be29_056d9d76c160(X,~,~)
%TP664F0162_7FC3_4453_BE29_056D9D76C160 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 07:22:44.
% 
% [Y] = tp664f0162_7fc3_4453_be29_056d9d76c160(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-3289.74882371951;-1498.84051842406;-977.81800685793;-2473.29678142316;-579.957172766836;-720.151909430027;-474.399996307268;-522.242249232241;-420.226401116416;-417.475292942104;-306.200888877846;-388.367090162115;-291.395845074123;-376.411586177093;-312.430500619519;-321.397730757885;-231.05260012331;-238.536185311164;-300.292838284052;-183.531087713643];
  x1_step1_gain = [0.000377693156370527;0.000577456077299579;0.000763241531595339;0.000656175798642298;0.00207889551360641;0.00132289146378878;0.0021050120182425;0.00188095402340317;0.00237297796658717;0.00258788322162195;0.00300539389135356;0.00292429008303262;0.00346972867854875;0.00273531031890042;0.00311974005335498;0.00346508516159047;0.00478122411926411;0.00448287798328093;0.0039349006097199;0.0052117714716069];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.2400617222607248635;-0.17644718363718070453;-1.1125715268561455762;-0.87493488409275388396;1.8565211605616791513];
  IW1_1 = [-0.75623416479767924248 -0.16204336539968325348 -2.1205595194398236814 0.85429581422389522238 -1.9110993063930510338 -0.066469045269775137164 1.1076673546434274709 -0.37327530191213531552 -0.98464205845315622145 1.5924001287341089661 -1.5308102935951837242 -1.7800220389179814973 0.1987858060256362247 -0.37614744076096184866 0.41338071757270472562 2.0990128318916361216 -1.0807380028854356446 -0.91198300930849951484 -0.060565364833220419749 0.38606579078257197146;3.1041449260941913479 -0.51202033829175652713 -0.57547618642782982779 1.7817633609423750851 0.62610882227945841016 1.1461530796736958404 -3.3475544627537243692 -1.2005366045330714275 1.5954004665271053209 -1.72387065979450127 -1.6030206765989891693 -0.12609762428251672373 0.12028969485060764244 -0.10167107398224480008 1.2303986495480625507 1.0680897382861069822 0.35207926100609332476 -0.32349593588260699972 -1.7312728595883009053 0.13165835321120145096;1.2414924412263490527 0.67208435026383506194 0.91222537478915954878 2.837070946025157081 0.34909271890329512944 0.55973480377229400951 1.0610065112193540759 -2.2091824323898179294 0.022907257508308177113 -1.0916217350037054068 0.58149603945681260786 0.77653130217294363824 -0.15602036673764338448 1.2445321828548407961 0.18272860616612263129 -2.0438076315928972448 0.35263277957388089456 -0.37477473036239811188 0.33733844603667523243 -0.76924273358978922932;3.5876576507860176157 2.5924453266453983957 -0.32654305264707378376 0.13201590058965487207 -0.976290799078566085 2.2882334475920997896 -1.1451871267972670854 -0.59239061491342537291 1.5325832224415731897 -0.15313498713235657034 -1.6647400587343867961 -0.70452672022286810805 -0.83587731583978452932 2.4285978192512649443 -2.8283502806686406927 -1.1282202288836800097 -1.4958411283769195155 0.2678326986385015096 -1.0517475427727513981 -0.43882255505667594475;-1.6738817372810017581 -1.6313425522612561291 0.5856673058249425301 -0.1622152345304170018 -0.73227918328831420602 0.25794497548617811855 -0.053458174477242971501 0.91158083304081405362 -0.024633191459228993836 2.5465702993821492939 -1.154605330471208946 1.5829551940363841833 -1.3187072404559732064 1.2345403317680332389 0.86428619689200192333 0.14790324499753074172 0.66225705746425966858 -0.59098865042623094368 -0.077457767585981024405 -1.6185274627194508046];
  
  % Layer 2
  b2 = -0.79568520436187317735;
  LW2_1 = [2.4796252243765937173 5.2778624564963241284 7.3017901444816430967 6.3690844067915080018 -4.4186301708130351074];
  
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