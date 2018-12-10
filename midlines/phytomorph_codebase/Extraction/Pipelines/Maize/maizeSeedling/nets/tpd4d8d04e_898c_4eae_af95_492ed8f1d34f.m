function [Y,Xf,Af] = tpd4d8d04e_898c_4eae_af95_492ed8f1d34f(X,~,~)
%TPD4D8D04E_898C_4EAE_AF95_492ED8F1D34F neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 14:04:45.
% 
% [Y] = tpd4d8d04e_898c_4eae_af95_492ed8f1d34f(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1425.12385949519;-2639.1224350771;-1040.58146532707;-1175.93373113978;-1366.539512037;-713.224898465731;-638.058802539169;-604.796304556692;-546.419507096578;-435.830146318436;-301.630323244972;-332.822721465589;-407.901071240808;-247.40747632391;-310.614407522926;-284.621358963187;-260.199856451159;-137.271225529371;-220.006532848223;-273.249410813717];
  x1_step1_gain = [0.000226370837385414;0.000563347828528141;0.000948114745185446;0.000858032283205083;0.000730242383199755;0.0015171639926653;0.00146770964020617;0.00151546500149801;0.00206494173697169;0.00202475070577022;0.00326439309746543;0.00229997229518741;0.00248425272472343;0.00412752991038807;0.00322794766513108;0.00343797851616053;0.00328303385537585;0.00723300236088577;0.0044572748674084;0.00407527899100744];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.8790817598161211333;2.7952276833968019432;-7.3845773586679506195;0.89867190115066930556;-1.4961470424226623877];
  IW1_1 = [0.39730574811770069665 -0.24128995019397403432 1.1839156154782641828 1.6299931405581664645 -1.2899275700037662329 -0.17113633143233014655 -0.26992570491321543624 0.15739529322629650099 1.106958466262636831 -0.89327913954874771374 -0.60779624360175554809 0.022896358129132335368 0.38446036187346455115 -0.86203329445768672734 1.3436745552028535045 -0.88201370752627139549 -1.0759495737720163877 -0.5368851759954137437 -1.5446605637308885584 0.63600031338578333351;1.8252189694227143324 -1.7165805053809064518 2.9314988243847914795 0.26936843616019184866 -1.0130895651210136688 1.5434666242521610791 -2.0274444374535063318 0.30876356614692390679 -0.48514297491592778266 1.824500354870603358 1.1755934790332076556 4.5368316528181171776 -1.9608554797059007946 2.1224899333111997279 0.88439394181297215081 -0.17842103400606340591 -0.7720413811630039147 0.14883602976901375015 0.60395542830107473709 1.214657313174752673;-5.1084158687232852358 11.061665431375702795 -0.37038079219361758065 -1.961757102398253183 5.9223106683189357113 1.6223966632296862311 1.1559365166042048223 1.869957724020487877 0.6356472035318840863 -3.0598833779014213974 1.7473075312088823807 2.3530936885987281393 -1.2910120059872016718 0.1767485464014196539 -0.92013649097942418198 -2.2694871940073046979 1.3195497075639943674 0.4673497704772466288 0.23757654054703913893 -0.16515077063071609276;-0.27555269400637061317 -2.7842392946587075642 1.9188430734937815814 0.62579539629669111989 -5.9797396376608711321 0.023804424272614758268 0.97194821258876618053 -0.58518852209228711114 1.9030388299028426768 -4.5455605790025446211 2.8020481767263567541 -0.59317267810812179807 0.42460906057751734988 -1.7496064657760981298 -2.9789722155066620246 1.4483522595672331246 -0.73258249559548394014 -0.24960136437175439728 -0.0010850540963457076804 0.59097310212882592229;2.0735935490299519657 2.0467623928101801312 -2.7499728421868860906 -0.74285586745071463532 2.9294279546211896559 0.90157376963711410855 -0.66014954798221725518 -1.3620822734912536589 -1.9037065746267813893 4.6491870027925799036 -2.6504111314550162781 -6.510978628404590296 2.131635242408438824 -0.69666843614951390329 1.8174529639765679701 -1.5337219880870320843 0.70990258876470269112 0.52993174649614593896 -0.75568284188881851282 0.46505926846638007133];
  
  % Layer 2
  b2 = -7.6254375780389027639;
  LW2_1 = [7.3319946296205724678 6.0233202825155345295 -13.307038487564886609 6.9972715807967329482 -8.0523395480963255721];
  
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
