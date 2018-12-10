function [Y,Xf,Af] = tpba9991dc_5dca_4456_895d_c6f02d1b1bba(X,~,~)
%TPBA9991DC_5DCA_4456_895D_C6F02D1B1BBA neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 12:45:43.
% 
% [Y] = tpba9991dc_5dca_4456_895d_c6f02d1b1bba(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1710.79601884287;-2767.31086728677;-1769.06834334395;-1424.63292253516;-1115.63541033651;-492.311800916002;-695.077686466816;-421.625702789062;-721.553885286555;-557.274532479575;-502.132109545491;-541.944235264823;-462.745727182158;-167.207345824944;-448.973304658431;-294.401776402714;-391.010406370086;-235.625647900174;-174.29482448455;-303.609020442431];
  x1_step1_gain = [0.000202107905561193;0.000554219398814734;0.000742185534904078;0.0010185341647551;0.00103332655312912;0.00147681597144009;0.00118221307527392;0.00265958726683864;0.0012523589192447;0.00165285753177999;0.00232860128237181;0.0018801318429768;0.00198822427969292;0.00479968082327771;0.00252047169225136;0.003631172061584;0.00272984263620491;0.00383651951608925;0.00502254586493077;0.00330771328744497];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [0.10215759404583080183;4.2093372105996786914;4.7962709805656773554;1.2725533731382150382;1.0188976779185441313];
  IW1_1 = [1.1486654246645162658 1.0171418291319787031 -0.38967743441817581651 0.91799940926775325245 0.14485778211660621517 0.55643582858717133277 -0.63920648146400094181 0.56939812315896920492 -0.86457580763453678596 1.8116940239406056357 1.2845004781814928219 0.084949911719935183863 -1.0867767930652483255 2.8642395949797521126 1.8498498741275737345 -0.60320940136912204466 0.98369920717557401968 -0.58848958539413520619 0.21529900896256670695 0.1232150469424269229;1.3323570087246325855 -1.3586044258876512902 -3.1658768557096479412 -6.0436236988102400858 -3.5614216486788414784 -4.006072907224097257 -0.29662874694707325185 -2.4505468848274976068 0.96791349744506405273 -1.5898958115208412512 -1.6279785271986486617 0.86890356486647024603 -3.0735856507438468199 -3.8171949670602578486 -2.3245390257306386417 3.0355046729842580788 -1.8779442365648357782 1.1035113942567780221 0.33242826838917594046 -5.4195899762208261663;4.5411155265905946976 -3.7948683738118513809 -1.141359624230766201 -0.99354201351957027732 0.4871890489527841428 5.2886912658716189739 0.78585259990625733195 3.8794794785553001759 2.2447066811860598712 0.6761169161445319542 4.3189292990790697502 -2.7086603169809047031 1.4756396371965778691 -4.6055624664231995524 2.2390057346338223532 -2.0179032154271387967 2.8300737370838486839 0.44063514141282900116 -0.57648228645193233532 -2.8809614748336125345;-3.6197551403621024413 -2.6309093258792701775 -1.5043900977511466266 1.4628881977485261867 -6.2265785868801666325 0.75306472894210607638 -2.4712559332045889171 -1.8224334674317201888 1.2701656720370593234 -1.320062103781947771 -0.45579973153592645563 0.24658161293373392908 -2.0042053786765627876 1.0821925425136871279 -6.5260505636128414153 7.5363972024652605697 -4.8393380378757315796 2.2528683448272190937 -0.053703323920149918791 -0.86438516199906589854;0.2393766009638479686 -0.61850848321809503982 3.6074655610434258612 -0.60309100997133058275 3.7320711396863379683 0.3222687864296418514 0.85352167353453956 0.41650127566030731296 -0.84032858761687867499 0.58453684107823145055 0.22129879896643009363 2.505932570492920064 -0.59802765142257396924 -2.460588638427651631 -3.0711157707170815279 -1.1984689697770141148 -1.8182383253605571838 0.71439616502911629237 0.56401831658697965111 -2.2579029365875258328];
  
  % Layer 2
  b2 = -1.3244009267524854412;
  LW2_1 = [-6.1012050720910542267 -6.2057799076851054565 -10.318458293827299954 -9.2997397356858932937 -8.5311078475687622813];
  
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
