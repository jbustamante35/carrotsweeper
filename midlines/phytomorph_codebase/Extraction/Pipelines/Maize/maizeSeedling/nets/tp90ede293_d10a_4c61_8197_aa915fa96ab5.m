function [Y,Xf,Af] = tp90ede293_d10a_4c61_8197_aa915fa96ab5(X,~,~)
%TP90EDE293_D10A_4C61_8197_AA915FA96AB5 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 20:31:49.
% 
% [Y] = tp90ede293_d10a_4c61_8197_aa915fa96ab5(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923585;-1028.71862143182;-451.710971542894;-692.538815519998;-441.043793232211;-783.667639824728;-450.026271289933;-656.555828069285;-380.160028036842;-183.589471391053;-268.743398235253;-277.584610398823;-242.574763141922;-313.061223281901;-302.321350914747;-296.419908591256;-275.821003918801;-162.786711803943;-176.062817024286];
  x1_step1_gain = [0.000192953812201724;0.000635170113819407;0.0010737117105794;0.000992967389065705;0.00140247416929523;0.00179426803546528;0.00152094087103248;0.00262475201653713;0.00199415494766509;0.00251075910029925;0.00454433874041166;0.00378662000738688;0.00302237737948937;0.00448879604138519;0.00366080733460372;0.00363057230533728;0.00381924645546747;0.00427537582274507;0.0059212055178134;0.00602422138644359];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.8028276639874951659;0.84768347228735874133;0.73454422168733179976;-0.79150569836606476937;-1.9809899686811154318];
  IW1_1 = [0.32734255240448406621 0.56491302045948343924 0.11588066630937061963 -0.10814964817548279596 0.72036358188066362818 0.10432259650497041625 -0.21575483412206414302 -0.41586537140802254386 0.017628709787708300716 0.064535576621952694953 -0.46304050275133579184 -0.71219485762644063609 0.44407269457107900745 0.5285414278173918845 0.23620697102386903343 -0.21599811292986947908 -0.40540700987387695564 -0.46671469987377706135 0.0021784538798495026168 0.3840474953315681228;-4.112759161449067058 0.8155029582531625465 1.2503903916265159957 1.861157107445888137 -1.0479612387622290193 -1.0936648280088090424 2.2053929873605015999 -0.27032019471083168494 1.661365707258124047 0.26877684815607583158 -2.2092744627784304257 -1.2864492291884941366 -2.8558062742837764247 1.1211080266790161009 0.28755167422018024226 -0.097634324181674114751 -0.32835257335291911929 1.1010993504842987445 0.54461824418125015512 0.65146355672976552498;-0.98128464200470555401 3.4158602652846883352 2.3158053991605775934 2.9343657505561666632 -2.2823676113918258679 -0.61054902094221386033 0.96606496755344861693 -0.29741771835330099716 3.1121598700267072424 -1.0171978770523277813 -1.2601769057820157993 -1.1517840253979987342 -0.4711962560912044129 0.17697637309448274734 0.84673713634887814994 -0.60623643083957057254 -0.51729491928934967504 0.048658011619867214126 0.81217804991138020654 0.078345592285197482307;-4.5408012132714485887 -2.4908247525790869759 -1.6575031008587368664 -1.6185135054970243296 2.3299458046402810041 1.1776794850402680925 -1.2463697586903816905 -0.24967370697477100916 -3.5088864191756954369 1.2099160426274504498 1.1526982937933054618 0.69117563628599287817 0.37353207296301471629 -0.64222595397227077996 -0.94360136701840957585 1.292889296706351665 0.45895847637837999811 1.2690560999030207068 -0.32984190940111340895 -0.43858727040980283185;-0.4094779634799234902 -0.009303701321383956932 -0.5708872131864671795 -0.037496101460639287561 -0.037758445837401073686 -0.48862515868376454886 0.31521564004351398935 -0.38728984379596886889 0.38353070548924139693 0.28277954209249628326 0.10979878179223820367 -0.14974953034110330408 0.54954186019323192358 -0.58956678222083991336 -0.022416943699461792244 0.038513809534663663581 -0.38009264748822357838 -0.037285035455518283909 -0.29657506391399168777 -0.36776161005420815453];
  
  % Layer 2
  b2 = -3.1768263073074507474;
  LW2_1 = [4.4156504686516564462 7.0179501116591369581 -5.4115512583242049161 -5.4978054545503196593 1.567487933905837405];
  
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
