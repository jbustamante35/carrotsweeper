function [Y,Xf,Af] = tp2424c703_e293_4183_8e68_d21fd4a11d37(X,~,~)
%TP2424C703_E293_4183_8E68_D21FD4A11D37 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 22:53:24.
% 
% [Y] = tp2424c703_e293_4183_8e68_d21fd4a11d37(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1710.79601884287;-841.368302509967;-925.675633657014;-1424.63292253516;-1115.63541033651;-861.953076142267;-996.664725988394;-421.625702789064;-721.553885286549;-557.274532479577;-502.132109545493;-521.810951623249;-462.745727182161;-167.207345824928;-448.973304658444;-294.40177640272;-391.010406370088;-235.625647900175;-174.294824484562;-303.60902044244];
  x1_step1_gain = [0.000202107905561193;0.000554219398814734;0.000742185534904078;0.0010185341647551;0.00103332655312913;0.00147681597144008;0.00118221307527392;0.00265958726683863;0.00125235891924471;0.00165285753177997;0.00232860128237178;0.0018801318429768;0.0019882242796929;0.00479968082327789;0.00252047169225127;0.00363117206158391;0.00272984263620488;0.00383651951608921;0.00502254586493066;0.00330771328744478];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-7.4868961270519855944;0.74507327709373716296;0.60616858663604300528;3.1911972987685226499;4.0269725390127790732];
  IW1_1 = [-2.4616397148758122349 -5.0065150241169700251 2.5959538975036751118 1.6722421969508538542 2.1439430158992718312 1.2016987211556113824 -2.0797853183311620917 -0.92322511051583078245 1.3868569594291999536 -1.0335270418302469064 -1.2051861605747453687 -1.0769253970810217513 3.7613411415192770093 3.1896161639928135578 5.630503029615788968 1.0882834578561533245 4.178986354722841412 -1.4180539640882470209 0.26731827395931911928 8.1521976280100201251;0.73655221938477744015 0.45162433041962302127 1.9075865586585565303 2.7398060275646436423 -1.1731279832554637466 -2.5076208006641111936 -0.23817136743021180556 1.4043208206787207715 4.5533247578751012696 -1.9260286881707819884 -2.3649139397123928319 1.1689935745189665717 4.408246257047331973 -2.8321224016377413157 3.4144198810589854176 -1.5685245278644102118 1.735032390712331507 -3.1300475961221634336 -3.1592720501686413215 0.024509433606589435539;5.0141888951916904205 4.1996241103100828695 0.33530054606424974928 1.5813993740585059466 4.0826265636921235114 -6.7505443621070693183 2.501946497181822604 6.8079516220126405912 -0.73978332574273830691 -2.3602474550768404171 6.2753253680455243213 5.8317159967204812432 2.3712790023984053889 -8.1205513622019402931 2.2121214282244130445 -2.2057155197909721345 1.8072343901544920541 0.15358534479150534602 -0.6048493064076125858 -2.7139382868633412649;-4.8163521420087338498 3.7713409824516306301 1.5749098163039891851 1.0682183650928331886 -3.7986679358098864867 3.1547123340742455611 -1.8542720109391555017 -3.8057112868036542253 -1.1697258981809934841 0.20473280742992919223 -4.8420095524659263475 -1.4522575446087953033 -3.1932029730905426845 0.34295652076615662507 -3.1754151385707451993 4.9489584072917374158 -4.3447517757081799417 0.54330143090336358025 -1.5957531341221942167 0.90916231816707016034;-1.6580218216288635791 2.7125934580369190918 1.5116093836329673739 1.8586753396143353978 -5.606628214906797858 -0.90141715031617652976 -1.2967328164374918309 -2.3802319940397844711 9.0637068110434952928 2.4720423887528193418 -1.6959354791981287036 1.7769020796671279872 -1.244925948610961397 0.39291583258777440513 -4.6334304776631798362 7.115287733042777063 -3.4386615443278905069 1.0193438919370017004 -1.1196529766125642436 -0.57870045245451406934];
  
  % Layer 2
  b2 = -6.4839071613246321846;
  LW2_1 = [11.041282579982501844 7.7736874722158075812 -12.099930254830830734 6.8627243813337788225 -14.850335341616647256];
  
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