function [Y,Xf,Af] = tp73d4321b_3bbc_4089_9a60_481fa29ef70f(X,~,~)
%TP73D4321B_3BBC_4089_9A60_481FA29EF70F neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 05:59:34.
% 
% [Y] = tp73d4321b_3bbc_4089_9a60_481fa29ef70f(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5310.24104946905;-3062.40716175129;-1831.68761111049;-751.778745457494;-800.631435004028;-561.234234396124;-567.818719457229;-566.708636588774;-460.813029179955;-365.327762675254;-365.703888775152;-497.927214156978;-445.050021503468;-352.492361124264;-290.108959264268;-350.309427866393;-268.266206223885;-161.912675664692;-200.676357162356;-192.300211130365];
  x1_step1_gain = [0.000187746732766858;0.000517587446146427;0.00076092913976903;0.00107683880161886;0.00113410115771159;0.00154637953830619;0.00191937835803698;0.00192268942981793;0.00211835684506625;0.00196958059572073;0.0026895539866255;0.00257377430529598;0.00284234016821655;0.00279908218986808;0.00409640929201516;0.00345584818973337;0.00422505290367575;0.00549698021899248;0.0047143261826188;0.00476878567270201];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [2.2641058226532124742;5.268305241919864379;-1.9580477716356192719;2.94176756897785463;-4.045540285753946641];
  IW1_1 = [-1.0381468651102128131 -3.3840556458753292546 -0.59421304337167724086 -0.40243499336402871469 1.9062102935489608768 -0.46866818164211127007 0.87074325424787657468 -1.8103370100761748684 0.48161757257777687569 -0.47084916517152958226 0.090032689340578411774 0.31809657612270553217 -0.19816218650582756444 0.60079741288646848485 -0.075540143558029695936 1.0803193435977589409 -0.34074283326730198107 -0.2912483911940477177 -0.18972692629006007725 0.051131105078423229326;-0.62731337115486307532 -6.4229460313117021641 -1.9209634083200279875 -0.63287016035670651259 1.9625758228466503308 1.584357386031676862 2.1860813422038138931 -2.5086585368031260046 0.39530154859043564741 -1.9282265503580071631 1.9094808490216077512 1.1509967567515073661 1.0148865333938874578 2.0042752462886248566 -0.47938896484859838676 0.85263538149174533665 1.749537949443315421 0.030790101585184879995 -2.6906893010799279864 -2.3113310648160476646;-7.3167517753837687522 4.2533614642036976505 1.9167075077061002109 -1.2949351025441422358 2.6406945856160137787 -5.1253120933400353465 3.5363177965120451418 0.98796087209874305923 0.3815957733897141213 1.9966708470180567492 -2.6934270712780725709 0.010866437613774327861 -0.70386578356334550399 -2.5635070499950680656 -0.35994250718124515265 1.4083651942189383544 -2.4664798435128179044 0.11213676430182882571 2.062187059105589082 1.9878429410525808585;-2.5032584249290703582 -7.9654745050864077527 1.1113925558397217408 -0.57771292803331186771 0.85617137779183305479 -0.39722195625037248767 1.2010681103932894231 -1.8551478323732506492 -0.042263117280211773275 -1.1916545841167467401 -0.294392837636193061 1.2968223903577396072 0.85095470599690803404 0.61003252038297828275 0.2386764015985584908 0.68530037399366472428 0.35663953592586145058 -0.42153796507611263067 -0.34418909178333889631 -1.0054513802692006674;0.82162155977411066576 0.92537869121578564258 -0.90971211272662999558 0.39037285663848969586 0.095140561589726391212 0.36624411534879230956 -0.30758119386224785918 0.063559714859755117589 -0.21681615832716047421 0.57894516023276598737 -0.31807329053654015416 -0.8594623359355841874 -0.30485707000420303459 -0.41541226158480432096 -0.96422022803887341702 -0.95819860841211657032 0.073993773117503908177 -0.060077581843738736556 0.31933879353553484837 -0.064769839220200992047];
  
  % Layer 2
  b2 = -2.4276905065714542431;
  LW2_1 = [-4.0360795930056667657 -8.8852104628730419478 11.757352530594582518 -4.2359314136686725405 2.2318696372148183471];
  
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