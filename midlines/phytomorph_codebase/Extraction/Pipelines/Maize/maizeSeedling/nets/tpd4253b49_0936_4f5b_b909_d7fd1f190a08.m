function [Y,Xf,Af] = tpd4253b49_0936_4f5b_b909_d7fd1f190a08(X,~,~)
%TPD4253B49_0936_4F5B_B909_D7FD1F190A08 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 05:49:41.
% 
% [Y] = tpd4253b49_0936_4f5b_b909_d7fd1f190a08(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5310.24104946905;-3062.40716175129;-1831.68761111048;-751.778745457494;-800.631435004028;-561.234234396123;-567.818719457229;-566.708636588774;-460.81302917995;-650.116847383029;-377.91382997879;-279.14169817331;-258.595524643432;-352.492361124266;-290.108959264271;-350.309427866412;-268.26620622389;-201.923452595307;-200.67635716237;-192.300211130338];
  x1_step1_gain = [0.000187746732766858;0.000517587446146427;0.000760929139769031;0.00107683880161886;0.00113410115771159;0.00154637953830619;0.00191937835803698;0.00192268942981792;0.00211835684506626;0.00196958059572074;0.00268955398662546;0.00257377430529601;0.00284234016821656;0.00279908218986812;0.00409640929201511;0.00345584818973328;0.00422505290367565;0.00549698021899238;0.00471432618261879;0.00476878567270241];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-0.60549594086107072499;-2.6990414965735047836;-3.674912673517722439;0.64327414012896777518;0.45913476075895581774];
  IW1_1 = [-6.4933610489759718831 2.9092915984327571977 -0.41810040314532254113 -1.041381008332119773 1.2074810076242807888 -3.3008566796837839341 2.7119495390886805275 -0.89874502561275004808 1.0933508096608186566 -1.4248137614624518665 1.7440702074611318384 0.48673925795957845297 0.50065989879856898082 -1.0043660520724100316 0.045314073631705503731 0.77304464077596035665 -1.2506476633002008292 0.49517325416594287768 0.084489642281184593187 1.2063878022355529751;0.23090780348532055632 2.6602252393947503073 0.20094163133812709132 0.29082692802766280105 -0.2592083651714282011 -1.1524731925453330916 -0.26111528227034769278 0.82383150344212019434 0.72442556229265697709 -0.92460008487982003889 1.0176890188248579605 -0.033172279721517074214 -0.42727015893659020529 -0.77825220030599107091 -0.1119480991475661491 -0.68904756203742734133 0.18147614248888579991 -0.16653018709389882934 0.91355500846345971677 0.077386385629612017278;0.104272393118422943 4.3176279712220804186 1.2554778557614527656 1.0959531770380046556 -0.31931177499764124006 -0.21186280205706739532 -0.030560016773392525979 0.95110511201068681508 -0.75477975958802823175 -1.1989626277428828693 1.2848710197673418598 0.71015521216050603304 -0.14003807451103039705 -2.1656402217216035311 1.0520708507310028157 -0.30546569971492992801 -1.3453541136109234966 0.53159595407591420102 1.8054297255926139432 1.5665740519923792018;-1.2510192556181172208 -3.4141475450007119363 1.349059095790362095 0.035541955332396041323 0.73619182968206631212 0.36339734390641753592 0.6955528858236028622 -0.43073953633123007734 -0.16155910676890192224 0.12996772646718737776 0.18278007361740514214 -0.54783098947312569571 -0.44393699335985598031 0.033219504271093498737 0.17596595612075055515 0.66150063851769169698 0.559165689936642063 0.71585602149304561337 0.08594752714526911197 -0.53897106534501448838;-2.8517624273789450307 -1.7828969859518697394 0.35329069312618316401 -0.26567090483883781804 0.92587081578510377344 -0.20842904338433565492 0.49240998726939999575 -0.39211530015660733817 0.47044519457318617617 -0.06637578385406078807 -0.43395550422184392403 -0.65699133522968133914 0.37640851736143621231 -0.013783524740438184847 -0.076041226623013449815 1.1257289594833552737 -0.36104282768723838792 1.2660068557425392388 -0.7348735998190526475 0.23816756557542093087];
  
  % Layer 2
  b2 = -1.8262565030090902596;
  LW2_1 = [7.7800600006287172761 3.5579860555063329564 6.5729892957823299326 -1.7851174484236009832 -3.8926550144400362541];
  
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