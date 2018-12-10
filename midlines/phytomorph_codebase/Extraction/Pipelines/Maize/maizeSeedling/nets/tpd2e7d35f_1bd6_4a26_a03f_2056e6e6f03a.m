function [Y,Xf,Af] = tpd2e7d35f_1bd6_4a26_a03f_2056e6e6f03a(X,~,~)
%TPD2E7D35F_1BD6_4A26_A03F_2056E6E6F03A neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 05-Sep-2017 06:32:43.
% 
% [Y] = tpd2e7d35f_1bd6_4a26_a03f_2056e6e6f03a(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079663;-862.937878923583;-1028.71862143182;-451.710971542898;-692.538815520008;-441.043793232209;-531.307871775734;-450.026271289932;-346.375264292285;-380.160028036841;-183.589471391047;-259.432124024266;-384.146123020662;-242.574763141878;-233.266353450978;-302.321350914855;-227.243600249852;-191.974128700824;-162.78671180391;-176.062817024308];
  x1_step1_gain = [0.000192953812201724;0.000635170113819408;0.0010737117105794;0.000992967389065707;0.00140247416929521;0.00179426803546529;0.00152094087103249;0.00262475201653713;0.00199415494766511;0.00251075910029924;0.00454433874041185;0.0037866200073873;0.00302237737948971;0.00448879604138534;0.00366080733460375;0.00363057230533655;0.00381924645546872;0.00427537582274614;0.00592120551781193;0.00602422138644332];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [0.42298987539757243459;0.32274332556208384926;-2.0310967095804808302;-3.7908860553255205872;-4.3400097562889197178];
  IW1_1 = [-3.5982935942346392366 0.4469117465621771923 0.72710612387827944225 0.23399205220585531095 -2.1606870539266336273 -0.16122383138711293693 0.10717295739871730842 -0.10959788633079951825 -0.54594221033123013775 0.24080271176138173805 -1.0840278465364290295 2.2179644423507109963 2.1303097425832864253 1.0478667393464522473 -1.3496329864275371513 -0.013740551981329261061 0.16906827992750372847 0.087592829513273878783 -0.81077251541484895991 0.42570054904744358826;0.11412543458847985989 -3.1129675246423036405 1.0458109791709366831 -0.22459971797570851826 -0.56208408387409625195 -0.17836019940874373035 0.57668055598654632821 0.64226108059666997541 -1.003665305057644952 -1.7549478073137945167 -0.56562718773137121708 -0.11763993353835708322 0.30481864880677472796 -0.62596354178908575339 -0.10502293654861744043 0.80351796415998366641 -0.019959508291517647632 0.28269382656145830568 0.2075623574214469913 0.49693598810837391611;1.0871326897748787577 -7.1348204552341947249 -2.8717554583120499245 -3.416951041429533209 1.8888857201175344169 1.112062688646795916 3.8624168970724941374 1.6676206332052605497 4.7469853492933982153 3.3290850842181929359 2.0201665210910224602 -1.3860072671105478914 0.16575613244088074905 1.212936824622116605 -1.9736829284086698788 -0.70461240864292284947 0.55299141204268364991 1.3258346227277169049 -0.48071425395500771982 0.073990912077425474713;1.1800929104321542518 -2.4789681489909916046 1.4156615097558780114 -0.21647423274312460606 -1.9682689012876637413 0.079265875742558500328 1.0175088238663247964 1.5476025255866292518 -1.3174862790478254482 -1.4069931666371702228 -1.0312318173652816888 -0.0061227396720530818941 1.7689808916526539306 0.17274333782953818495 -0.97667918916262952411 0.81351763139296007221 -0.15424763297552882002 0.86837161417478703207 0.28281290008741860609 0.55443086213689685149;-0.27059418013951347382 0.44131099663557910029 4.7402807757756608709 0.91827380487387488817 -6.1459976836322134375 -2.8680620370737663904 -1.0308255204489447099 -1.2329427235948207109 -6.8779104274197093716 -4.3079687248774751751 -6.0361636054299259158 3.3780012298420909467 6.7784543057824144796 2.1150316747296931119 -1.1575869693852500752 -0.0068398398136126879743 -0.32938105073317869298 -1.8390976554941256271 -0.27423407523464388147 1.2352121513253924157];
  
  % Layer 2
  b2 = -5.0664322984057514887;
  LW2_1 = [2.4313772006899601053 -4.0548950284868574556 12.611466254218210992 5.0991440795041826917 12.132924866917338846];
  
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