function [Y,Xf,Af] = tp2af94e9f_fe4d_4642_a76e_005db9226b15(X,~,~)
%TP2AF94E9F_FE4D_4642_A76E_005DB9226B15 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 13:47:13.
% 
% [Y] = tp2af94e9f_fe4d_4642_a76e_005db9226b15(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1425.12385949519;-911.082072190464;-1040.58146532707;-1154.98089660502;-1366.53951203701;-713.224898465724;-638.058802539176;-714.930642633132;-546.419507096587;-435.830146318432;-311.040988170551;-332.822721465577;-407.901071240821;-237.143827584618;-310.614407522929;-284.6213589632;-348.992764796145;-139.239122485351;-228.69812627522;-217.514535510738];
  x1_step1_gain = [0.000226370837385414;0.000563347828528141;0.000948114745185445;0.000858032283205085;0.000730242383199753;0.00151716399266531;0.00146770964020616;0.00151546500149802;0.00206494173697165;0.0020247507057702;0.00326439309746566;0.0022999722951874;0.00248425272472334;0.00412752991038843;0.00322794766513109;0.00343797851616046;0.00328303385537576;0.00723300236088659;0.00445727486740867;0.00407527899100751];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.3479640914517538874;-7.8109734286098966294;-3.4398900902626912135;1.9596964141268784498;0.95316441082443892352];
  IW1_1 = [-0.38496189893516707947 -1.336130547410768088 0.39867809290561234681 1.7832037775211559971 3.4059731930360590546 -0.11897265320926214693 0.60888722326278610364 0.27916721619237372387 -0.47844884479462723714 0.90041235736002445122 0.18863360659990768831 -1.7324183286000798265 -1.4579575318350896929 -0.81019059964592932133 0.91495686400469489197 -2.0197315234710071508 1.1382583931285215062 -0.26252694605409959694 -0.30617745241424904057 2.1930073977062858326;-5.8948386202605469464 -11.326557649712594156 -0.88225217193519078762 1.2695349018160619892 5.9926613373617216851 1.7202489170580919087 1.2993276942055707046 -2.2669962001080210534 0.62722196026698340265 -3.3281662282132793429 -1.7050911974895173184 1.9488337887403088811 -0.090385539270900427944 -0.24384703957871289859 -0.69236244075720765334 -1.6988479950505370741 -2.1883746964733279405 0.070582918158375862472 -0.029774288043216440292 0.79624984536722509532;0.066979936242048732087 -2.8995841879604205005 -3.6452992885604418305 1.9228694558656949098 2.6046960142309916186 2.2511071173694872094 -1.7012256290349316856 2.3326740479823455665 -2.0300640084465975121 2.5378109091032556854 2.0470124197342238759 -4.9021070638318846591 2.0423977040427589991 0.62652245437804965356 3.7642219369118836703 -3.6471535138134902887 -1.4896111086353531761 -1.2526761563278039358 0.34975793635233443091 -2.260405077432962262;-0.44043156820620876779 3.4810122531495624187 0.4307723210908118161 -0.14232400561245039938 -5.2250282914188526107 1.5616319687720587872 1.0969235712494398438 1.7030474752301782981 2.7748179046278473869 -3.3023521670318425336 -2.1566976173536933636 -0.21978090344495659902 -0.56825744536576816568 -0.24954033227905966785 -0.539371303666557389 -1.8124821473186372156 0.039064996088580127009 -0.25793429708970411207 1.3742319275895398256 -0.58549543474789533803;-0.66246808845109272923 1.5749222481913736083 2.7257044786968034344 -1.7554126980491524535 -3.0135180175897269272 1.5934594588971033868 -0.087628250153274481193 -0.40433348924734657803 2.0922451217080522134 -1.3376883248750450761 -2.4512377470903379262 2.5507832502806366826 -0.55743626806837298648 -1.1920658477470300607 -0.40158051585906706649 0.038626073025936159822 -0.56869348862364310637 0.47976492813241333746 -1.0117591108824566248 1.4429893870820038604];
  
  % Layer 2
  b2 = -8.0954233511020010639;
  LW2_1 = [-8.9588229322489212336 -13.326561124384987522 -7.7203547461683070807 6.9704483228722438071 6.1033687402140790468];
  
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
