function [Y,Xf,Af] = tp2696aff7_a109_4f9a_99b2_f45b218c0f76(X,~,~)
%TP2696AFF7_A109_4F9A_99B2_F45B218C0F76 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 18:06:55.
% 
% [Y] = tp2696aff7_a109_4f9a_99b2_f45b218c0f76(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-839.604414481266;-7316.6928526506;-873.097427232645;-1646.05375716083;-423.333970720151;-346.884767278008;-410.542636124639;-450.689400777709;-342.30806526671;-280.547165434103;-739.082929249925;-246.980329538456;-167.382252836757;-227.785437026333;-141.910346711091;-156.623146275721;-201.347350487964;-147.229750402561;-106.365322334024;-185.847378023487];
  x1_step1_gain = [0.00043525835929838;0.000254266760874674;0.00034735300303884;0.000777513986910232;0.00172505617054034;0.00206247464333953;0.00302249934878054;0.00323098252238321;0.00220614635584281;0.0024281208277307;0.00231897377235314;0.00407495867726348;0.00341201206501292;0.00551378415971723;0.00648970220781294;0.00620026189110773;0.00633602455501666;0.00774864116842937;0.00243421857993183;0.00313801600663448];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.2908035673347832972;-0.094377254179652911303;-0.23605315469810925211;1.3987589032192535576;-1.4738486108249275741];
  IW1_1 = [1.430411846365898354 1.3233481131133946196 -0.65311609378688872951 0.97678506635702933014 -0.91325409149016234522 0.13641864137006620705 1.2257076276035412388 -0.039716075955611690462 0.31034526613624158209 -0.21609488645474264623 -0.053283022958475950714 -0.28723736950802947199 -0.27322893864748204473 -0.23243874774071915379 0.14663177693968038051 -0.29549990474604881374 -0.1000532711229176891 0.075528043828272259375 -0.19058030413969070938 -0.2382050676766760966;-2.1524975021041781176 -2.6044789915038157879 1.2522343363924759441 -3.2425787987889975739 0.048619111614162463519 0.33138988739650060422 -0.80809600301367856545 0.23490267071700224322 -0.23215857907628320578 -0.30912342338071724379 0.98465013333090156689 0.10837376167741899491 0.023768525568174103668 0.68050604523541990964 -0.035760568512852593248 -0.45488306439326597896 0.70107714340264348785 0.31422065446263813815 -0.63984133695385125851 -1.0165599860220633222;-2.9699208419288414262 -1.4669670258442393251 -0.48043410702969602255 -1.4792096351472858728 -0.54833456785470269246 0.61173971674385463349 -0.29447325743854674274 -0.5832178703258533714 0.4656345613890928159 0.33645732142573853141 0.159285095745903299 0.077275260441227996977 -0.38913350571571603265 -0.31311469620183718376 0.062814639817365125474 0.07782479162384359217 -0.3215163210500211477 0.18289106378081773641 0.58103456304287615097 0.40961843529876174053;4.5532659301745681546 1.2638950035566312025 0.69078332010227572191 2.3310250588623588897 0.76302439363621776813 -0.41127365924305719913 -0.016216797108310337938 1.2102389461768690015 -0.71022578909242650713 -0.73470306246944938611 0.32493623713016622778 0.74089130217657339639 0.22579145461520241933 1.0819045403994811405 0.14921354101613826248 -0.97719415214051508212 0.12022997711054302095 -0.57253561419427445589 -0.71811479886906293135 -0.64606148991410983218;-0.5414192692823348585 -0.3835181852544104486 -0.36802366429991689722 0.14184883267544687313 0.13304381141782992137 0.22163653307943867454 0.54190077467670294542 -0.28235863881899370131 -0.42934305385702115077 0.45553650475185752944 0.35394996591280780684 0.063280661966308743294 0.47081230538848745093 0.31491019189310254234 -0.12858150723837691309 0.11528074097463730796 -0.28395138521597695913 -0.43633360944499927081 0.21772773898253719449 -0.49223616081457188809];
  
  % Layer 2
  b2 = -1.7542926538309089324;
  LW2_1 = [3.2944388069527725804 -5.2283295606443678238 -2.2559730724578073158 4.5812567446744933974 1.065026078339946336];
  
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