function [Y,Xf,Af] = tp7a6ae2fd_8880_4927_99f2_913d6a0e921e(X,~,~)
%TP7A6AE2FD_8880_4927_99F2_913D6A0E921E neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 12:05:41.
% 
% [Y] = tp7a6ae2fd_8880_4927_99f2_913d6a0e921e(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-8184.90793418462;-2767.31086728677;-925.675633657013;-538.973276654921;-1115.63541033651;-861.953076142263;-996.664725988397;-421.625702789076;-721.55388528655;-557.274532479568;-502.13210954549;-521.810951623251;-462.74572718216;-167.207345824921;-344.528961662822;-294.401776402704;-391.010406370092;-285.680184537944;-174.29482448456;-303.609020442438];
  x1_step1_gain = [0.000202107905561192;0.000554219398814733;0.000742185534904078;0.0010185341647551;0.00103332655312913;0.00147681597144009;0.00118221307527392;0.00265958726683856;0.00125235891924471;0.00165285753177999;0.0023286012823718;0.00188013184297679;0.00198822427969292;0.00479968082327793;0.00252047169225133;0.00363117206158392;0.00272984263620492;0.00383651951608922;0.00502254586493046;0.00330771328744481];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [2.7021844758276003517;-6.5536557217790667451;0.21932048917754576323;2.0100194168292868646;3.1811603347773904638];
  IW1_1 = [0.49171150074943087427 0.79320469564106432792 -1.3174186464919788175 -0.58522667815800866187 1.1232941393283659703 0.011872487848801776406 -0.40886114130394629607 0.85922788125576599061 0.77020405382541179407 -0.17811766469954171122 -0.9727667039122778192 0.27621038066537356714 0.027655304269224228325 0.63246178062428881628 -0.26040991265830049084 -1.5104415174176253345 0.22423130472622326836 0.87253106397156554586 2.6180539180211974681 -0.41930252334648476786;-0.63805219673297908223 4.8438184071476442227 -2.4721202096301722051 -0.09328533473277747512 1.2812697016475054212 3.199791226302308278 1.4009243027506488311 -1.972687104910289424 -1.3945015798788571448 -2.6013383618097858374 -4.8732348130621563342 -3.7386019945705859513 1.2385318291185107853 -1.6930092093797495245 -1.3478679451133965816 -6.7565284878970279792 5.443867238225655214 0.51938224104851726803 -0.82478808577504636457 -0.21619472256498126206;1.727033422028490639 4.162322889224291167 2.3168971786100951427 -0.50908483064367215931 2.0240592039647453682 -3.0330038693547090034 2.1901352543196721534 -2.4898678495394719334 -2.3260413670207116255 -0.095564964764970977695 -0.093732082910383276997 0.77745906477585535743 1.2766631003231978525 0.75996001270667090655 -1.0595041766302408792 -0.72534100128317613443 0.70674732966756725894 0.26654370791772646898 0.99228386592034378832 -1.9193939668897641582;0.74494077135282310564 -1.9971797809675346791 0.84345929919513062956 -2.6649898669249529881 -6.2297979854364857744 0.32102473065018932719 4.0321798752563848822 -1.0936202377122901108 1.5130900234918143177 -2.1972007491820786385 -0.71129017032709507262 1.7165913248074169406 -1.7845351649249956427 -1.3049040226187886837 4.7267487257561322522 0.89565285549014450606 -4.4488189248927838548 -1.4246705481457297982 0.4408280352985300965 -1.5080214873621373695;-5.8823721209114037123 1.2170621419356828508 2.0112148111643324278 -2.1689936361310722113 -2.0852359055443701052 -1.1925547157062637016 -0.59515629109874002012 1.7134431061720936107 3.6561195040722318517 0.58624587474675304843 1.714430240915269632 3.4250498201532719733 1.3346267450784636832 -1.6097271637362853269 -2.6533433369865169205 -2.0352270453481544621 2.0564273005258946192 0.98771684118321934065 -1.5538968811097495504 0.44394639035225219592];
  
  % Layer 2
  b2 = -4.5532163696813467979;
  LW2_1 = [-4.9430934062838840504 -10.715915600966628318 -8.2951114684433111535 9.355356065182473202 -10.434804659083436462];
  
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
