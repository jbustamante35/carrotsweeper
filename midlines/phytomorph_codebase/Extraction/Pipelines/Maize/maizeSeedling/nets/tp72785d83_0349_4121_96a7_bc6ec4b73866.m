function [Y,Xf,Af] = tp72785d83_0349_4121_96a7_bc6ec4b73866(X,~,~)
%TP72785D83_0349_4121_96A7_BC6EC4B73866 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 05-Sep-2017 06:59:38.
% 
% [Y] = tp72785d83_0349_4121_96a7_bc6ec4b73866(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-3289.74882371951;-1498.84051842406;-977.81800685793;-2473.29678142316;-579.957172766837;-720.151909430019;-474.399996307269;-522.242249232239;-420.226401116413;-417.475292942101;-359.26928651379;-388.367090162114;-291.395845074129;-376.41158617709;-312.430500619524;-321.39773075789;-231.052600123337;-238.536185311156;-300.292838284029;-200.215592449859];
  x1_step1_gain = [0.000377693156370527;0.000577456077299579;0.000763241531595338;0.000656175798642298;0.00207889551360641;0.00132289146378878;0.0021050120182425;0.00188095402340317;0.00237297796658718;0.00258788322162196;0.00300539389135355;0.00292429008303262;0.00346972867854873;0.00273531031890042;0.00311974005335496;0.00346508516159045;0.00478122411926398;0.00448287798328077;0.00393490060972029;0.00521177147160692];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [3.0857875356733819316;-1.267643728839229933;1.0620116320645718311;0.2927804759639833776;1.5068621780281350464];
  IW1_1 = [-0.56775712045154180796 0.19777643714416170972 1.5049932608056788563 -2.1600911516161196424 -2.8452243014199463111 0.94225084357583588357 -2.4161915820636674646 1.9040924812375097375 1.3874570316048244756 0.44098305904053614057 -1.1723723148220943635 2.1668601011753767516 -2.3591710098418197106 -0.23765752660201394653 0.39668898116826967204 0.89655100610283666729 0.96923472919336450815 1.0340914925555126924 1.4889007434267724772 -1.8010468744434151578;2.6182094144207113118 0.90222062054473550763 0.27551021067950137722 3.4970468762950455144 -0.99921753292981763117 0.1984388530535641515 -0.75739429863061569215 -1.2732937696588380039 -0.5270920946747100766 -1.338482161058640596 -0.9774324167708402511 1.760667184118436257 0.91098701900419132294 1.1493097781781982381 0.2989277165458046337 -1.1009542336482605052 0.67295481002439505591 0.036711236461509449969 -1.0980364644841884569 0.31600129054379649807;-4.0761111023584515323 -3.3510498612036507815 1.3107449590394590899 -0.12940457441462174804 0.95199282477381275136 -2.2462051411056878081 1.1367758335019766402 0.95526053209572914238 -1.6869801442145033743 0.62270146618204202937 -0.43741160253199856678 0.70246805650149901634 0.6977145925738128529 -1.5185155015264397704 2.0255460205736053325 2.3136171981641453499 0.7155589253263370253 -0.56407459778140411899 0.77975142682506504155 -0.91254509865383592881;-1.7488645995602676297 -0.68612449180596690734 -0.64011013791549675744 -3.2086108566451487789 -1.9326617667074659224 -0.16589273056827094632 0.71106699861537281393 1.0113091122233461849 -2.760870235166523301 -0.3613697724537477951 -0.25076151459007423039 1.1793411101057820201 1.9030522130943008374 -0.59009712231539135363 -0.78952857606510018496 -1.8808586648307648037 -1.3998741120325244314 1.1141944006256701005 0.94700686385990839877 1.0610063476812123806;-0.6185245160715567847 -0.579351101595618867 -0.26349858947904097883 0.12106431466868965152 2.1662741339683595498 0.6918889950796084376 0.11925880038836261798 1.6657660684069874879 0.32713390567244321483 0.42930960302109755222 0.7347908490501329215 -0.87136221930529877966 0.79922308648380557816 -0.030522436894403183472 -1.6379902748540411395 -0.34123449734687794077 -0.27889738039041689355 1.3493913545582889579 -1.1072698732549439704 1.6340950577788457831];
  
  % Layer 2
  b2 = 0.48794541181754153003;
  LW2_1 = [-4.303202895233741998 6.2301202915005617911 -5.365867793365923788 -4.961257669154655936 -3.6748996163287608319];
  
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