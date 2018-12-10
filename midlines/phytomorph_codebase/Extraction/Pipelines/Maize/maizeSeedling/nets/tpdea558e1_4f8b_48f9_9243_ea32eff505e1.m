function [Y,Xf,Af] = tpdea558e1_4f8b_48f9_9243_ea32eff505e1(X,~,~)
%TPDEA558E1_4F8B_48F9_9243_EA32EFF505E1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 16:00:16.
% 
% [Y] = tpdea558e1_4f8b_48f9_9243_ea32eff505e1(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5973.6023996224;-2528.40271623793;-1416.92748927618;-923.044887924526;-742.627689023617;-590.342170311173;-636.877630731671;-474.328050496929;-408.756773605243;-437.689914113497;-597.158358369105;-384.721680725394;-357.743301501269;-305.701917871576;-262.500890231163;-229.281637810815;-306.92449968068;-299.148398756998;-279.306010526068;-382.085749121026];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296744;0.000975455655635302;0.00136867839470428;0.00177765334299387;0.00147664952758212;0.00196526984478519;0.00211681540578921;0.0019139529940499;0.00213803411115013;0.00313745816598823;0.00239499394158855;0.00262802245354375;0.00340299377226972;0.00307580218448474;0.00420342200465309;0.0035393752665976;0.00305696135696249;0.00320422591704709];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-2.4857653567813442486;3.3021265325014321768;0.52450109104288433759;0.70105880738506898275;3.1605692414796440204];
  IW1_1 = [2.9882589110257922194 -0.34052727530575671588 -0.68443680148715746459 -0.95719255143393799923 0.56049917162474227617 1.5494131551872725172 0.27572525396719260726 3.7826786188638146236 -2.1504546780755044999 0.99775330068713319331 -3.3091722041004332233 1.7260271659648864784 4.2535085901430127464 4.3120888046348060385 0.79609781747351804349 1.7285992780469221497 0.8986115793585766065 0.16671506291883783635 -2.912298687045446588 3.7405078039225414521;2.434907356826566005 -8.2495654709208281474 2.7205684062789354805 -5.6235037547051183182 3.7010471019285344774 0.44191215084794349677 5.3219187493562554181 3.2148933257230085658 1.9549983895686307012 5.8370287146535648759 -1.1953153240319660533 1.2572207330580316142 0.95224660418884765622 1.4451524007493019575 2.0007578219485631088 -0.37299666240539763873 -0.20243386355532919096 2.3525836633316106195 2.8560800540100208522 -2.5037252199545614317;-4.1744202103014877991 -0.69051531006837596482 1.1116493786070051897 0.26024348320089824016 1.0308112135352514027 -0.8772956807472631624 -1.3171934109351020048 -0.92195421260694510046 1.2143208615423097818 -2.339977082902138239 -2.511107197859458573 1.2531831521639695826 -2.2435845534193239459 -0.18048047468468605969 -1.7740550505817545002 0.44081741369107518835 -2.1223107666220801804 0.73412584227119048741 -0.91619534784309608622 -2.1325782309943339676;-3.6552860008524388213 7.090710202690728714 0.27115247960660077142 0.58490034009669811255 0.46197817723309364535 -0.37263272425861115877 2.1529064391216281571 1.5414578561889944464 -3.4455158092291955363 1.1347190467431982341 -1.1263763431997442233 -0.96942703213983627553 -6.7163115965859923051 -2.6363536755339693762 -3.1330623570686264578 0.46707767076736295131 0.028329236649972066447 0.38620365582998811282 0.38293680104004002018 7.698407069355337029;-2.9593378792665130383 -7.8555903871706709296 3.1958067946965273975 -0.38919289968529668933 2.0558268608456775617 0.053999186168393303875 2.66692364089712175 1.529066532307381765 0.54698665679809010953 6.0414311186258871089 -0.8767086518317421806 -2.9296037525983997796 8.0776485499083641173 1.2666087951294224823 -1.3568883072331525685 -1.341696599611818419 -0.18994950988673342529 -1.1452034133967923246 0.054615402175969662846 5.2919850305371323174];
  
  % Layer 2
  b2 = -1.830085361767641805;
  LW2_1 = [-6.4597204541710171455 -8.9107538119702649482 3.2599209936536066579 8.3309398845183650906 -12.188721787835451238];
  
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
