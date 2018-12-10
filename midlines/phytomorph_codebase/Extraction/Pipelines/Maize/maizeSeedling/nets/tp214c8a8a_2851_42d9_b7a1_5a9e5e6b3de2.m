function [Y,Xf,Af] = tp214c8a8a_2851_42d9_b7a1_5a9e5e6b3de2(X,~,~)
%TP214C8A8A_2851_42D9_B7A1_5A9E5E6B3DE2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 08-Sep-2017 12:07:47.
% 
% [Y] = tp214c8a8a_2851_42d9_b7a1_5a9e5e6b3de2(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020144;-780.471896742374;-460.20555521635;-234.310604878461;-540.291178237134;-572.142715533789;-620.900292310244;-263.800490630846;-210.820779921779;-258.183122033526;-286.590370604324;-169.741999038739;-193.663740648594;-194.470619647338;-373.10470912035;-213.719449945769;-260.914631799833;-135.174215342508];
  x1_step1_gain = [0.000342540442445779;0.00027748074658588;0.000806013241167018;0.00156677380853905;0.00168393926281428;0.00269046527889069;0.0019808499204541;0.00236361816544284;0.00245784730269615;0.00195235680232595;0.00282106579415898;0.00424658079482957;0.00317261769260759;0.00555458687461896;0.00429724736478089;0.00484270636687461;0.0033393160809906;0.00598913897545074;0.00384783002604038;0.00642488839732234];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-3.4542853015163208497;1.9522500014678965119;3.8965635938815843531;-2.9827973638218181485;2.4170829765286505619];
  IW1_1 = [-1.6111354306220881938 -4.5120944424068829548 0.025748714757205713732 0.75707347708728811586 -1.7560877049521486093 -1.3984643409418404669 -0.14733904029884145603 2.938977178161480186 1.7647742932608971866 -4.6845215975808329745 -0.59927933953193690542 -1.6259407034392863167 0.92169010956640617938 -0.3533570670013272963 1.2425145766675447234 0.85610011657603901636 -0.13897507235357356437 -0.43734507500265218383 -2.233326502497535504 -0.31307246576953534634;0.85087225038676272959 -7.4457698405510708284 2.9990679493014580359 1.0965269845525564829 0.16693339667157475192 -0.97078139936948071309 1.7167929794928793452 1.1151054245868157189 -2.3125332484222287377 1.0967132909437546928 -1.7486225338140635177 1.3104942551569913611 -1.6076502003548200914 -1.0543226677836856631 -0.015749922785452327711 -0.22930659404217315633 0.33859008651127170397 0.072072955578903211138 -0.87765164827483621579 0.43872041871456246565;4.0833230933298860776 7.1551673734505687463 -0.53612064485072463693 1.1468867342340425264 5.5361844458087805521 1.7536145756752896663 0.051254009923396126602 -4.1369006140851816156 -1.4887445993302272917 0.94155397633133286472 1.9863580996953249524 -0.33850142241205777305 -1.9208956540276072023 -0.43881446266124346511 -0.37510069726894074149 -1.0607798936368708631 2.7319872390443884491 -0.70270296530579412497 0.45189512171649726113 1.275026590006796301;2.9407132433984424225 4.5963323514591829522 1.6913341868870956919 -0.7361290349106243136 1.1193136810370785117 -2.7138183256982006952 -1.4025511137393067429 2.174307396909354928 -1.4688762087490199182 -1.3053787628832660417 2.0807401848183291193 -0.29048511740025012573 2.3866512655007867671 -1.4356163731332238775 -0.88063253392515461471 -0.99339436968221028756 2.6910297408378243134 0.19401802838179069144 0.54416031716335766522 -0.78134706564213762814;0.032691554721118690208 -1.2045937452057215999 -0.95034120844265068673 0.35631427040562957398 -1.610836679821822548 -0.96415844890944257806 -0.25694489103697076837 -1.1924235895992243428 0.52338560208095263881 1.2003110949051343237 -1.0348963488356592855 -0.061252501067814317526 0.30749172514713862814 1.8589723728542597581 1.217070739056739237 -1.54868255747284711 -0.88621161502143830369 -0.79594561715198697271 1.4982139781805268175 0.088749777344230659781];
  
  % Layer 2
  b2 = -3.3616399594686878771;
  LW2_1 = [-7.9391264536855761236 -7.9058960193967120489 10.583563084617990668 8.8931609934647699589 -4.4083171527008433443];
  
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