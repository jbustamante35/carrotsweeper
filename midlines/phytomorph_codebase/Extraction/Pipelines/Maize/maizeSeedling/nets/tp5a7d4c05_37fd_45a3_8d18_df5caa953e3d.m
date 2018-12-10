function [Y,Xf,Af] = tp5a7d4c05_37fd_45a3_8d18_df5caa953e3d(X,~,~)
%TP5A7D4C05_37FD_45A3_8D18_DF5CAA953E3D neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 05-Sep-2017 01:20:05.
% 
% [Y] = tp5a7d4c05_37fd_45a3_8d18_df5caa953e3d(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-2001.76879548539;-2528.40271623793;-1416.92748927618;-1127.27896682562;-718.635970677915;-534.736579087551;-636.87763073167;-474.328050496919;-408.756773605265;-437.689914113484;-597.158358369105;-252.737018064957;-357.743301501239;-305.701917871529;-262.500890231165;-420.955237007915;-168.878309979169;-299.148398757012;-374.938438931838;-382.085749121009];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296745;0.000975455655635303;0.00136867839470428;0.00177765334299387;0.00147664952758212;0.00196526984478523;0.00211681540578909;0.00191395299404996;0.00213803411115012;0.00313745816598816;0.00239499394158852;0.00262802245354397;0.0034029937722698;0.00307580218448479;0.00420342200465308;0.00353937526659751;0.00305696135696243;0.00320422591704723];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-2.2347113907201689997;0.99621427123807149151;1.3340211961254844208;1.6477382211313573723;1.434951299697188043];
  IW1_1 = [-2.0458242479783446122 -1.6844431087694955451 -0.35280483508487497923 -0.41555462526129904877 -1.222210264446675243 1.0572368068567670196 -0.098898582640576992997 2.3900295087031597774 -1.3336720798556618295 1.8482981354589838041 -1.4071248876822530782 -0.56596588127249214395 3.2652198696149419277 2.1146214595918904067 -0.70174892958413970234 0.80422994590067187204 -0.2375370701596515588 -0.97948234258000510533 1.2193653508701156607 3.581173885971909332;1.7714846083319091274 -0.17471056994999398237 -0.36786688728463085019 -0.64016295728203131965 0.90133008475139941851 -0.71366363587406844715 -0.8334310016406517585 1.4226477593335240535 -1.8277287238925514767 0.75287298626476240049 3.225592974546499736 1.2281401256975037484 0.96312075920170170384 -1.2494214108903216598 -0.54933672645947362323 0.87466104926251520002 0.60802999029254811703 -0.28674566236959830867 -0.39920907660598881073 0.085038784376497866813;3.1850829102753990263 1.3768799792765573642 2.1392564984633879455 -1.5226432004351504457 -1.3279530496368641135 -0.81128353009335696377 -2.3148453383564620012 2.0676898516633301917 -3.1712080430897686689 0.99632721368813115248 -1.4051631412368372764 3.2808994474232742178 -11.442387204314977467 -3.9761452004169854924 -1.2160390133204006524 3.3052973694704395058 -0.9674782174859283268 -1.3465184428988734577 -0.4225410859319156387 2.8273128666355784588;-0.8159837612815138419 -2.9441160964545223244 2.3191562892263659634 2.2485573610576041759 -1.5111390133525748336 0.50799888118706881812 1.6144567795354891615 -0.5469789939552925917 1.1297431742689771905 1.6443858277264498113 -2.7310697801070094037 -0.41562953204433761289 1.0439418137148408494 1.7655042921364485498 2.1899009951293506582 -0.1153883372307998173 -1.0377496726242712022 2.1177491251294529029 -0.47113842885227957824 -0.87804616357879383681;0.71384662432267542709 -5.9692442991540097808 0.35080435827962369588 1.9254637717571199129 0.0080276188386336384817 -2.281009079067798595 2.6208924790915877345 -0.52990198470463301117 1.88171742620929483 4.6731273790828256054 -3.6773072262252206599 4.4859694607072579942 -0.64035520686552027669 0.46365213554694545417 0.19995240048803364363 -1.293917033019843954 1.2442926496550028581 -0.1734583227104621872 -1.5378989775602072143 4.2134722308550260905];
  
  % Layer 2
  b2 = -6.693261686430168389;
  LW2_1 = [-2.8324628213219793516 -5.7444087166153607882 7.4671103888500667267 -7.0389534727378304524 -4.676801760139260189];
  
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
