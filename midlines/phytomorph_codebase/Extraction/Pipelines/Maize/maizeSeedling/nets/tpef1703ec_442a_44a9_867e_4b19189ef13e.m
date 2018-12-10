function [Y,Xf,Af] = tpef1703ec_442a_44a9_867e_4b19189ef13e(X,~,~)
%TPEF1703EC_442A_44A9_867E_4B19189EF13E neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 23:52:45.
% 
% [Y] = tpef1703ec_442a_44a9_867e_4b19189ef13e(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5973.6023996224;-1846.28602050405;-1416.92748927618;-923.044887924525;-742.627689023617;-534.736579087549;-717.539895324647;-474.328050496929;-408.756773605247;-607.267828431681;-597.158358369104;-384.721680725411;-357.74330150129;-455.326587536788;-262.500890231169;-420.955237007927;-168.878309979179;-299.148398757001;-279.306010526063;-382.085749121053];
  x1_step1_gain = [0.000250772026915415;0.000457175383291265;0.000795694570296743;0.000975455655635303;0.00136867839470428;0.00177765334299387;0.00147664952758211;0.00196526984478519;0.0021168154057892;0.00191395299404994;0.00213803411115012;0.00313745816598811;0.00239499394158855;0.00262802245354363;0.00340299377226968;0.00307580218448474;0.00420342200465299;0.00353937526659754;0.00305696135696251;0.00320422591704701];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.7949126789287643824;-2.2591621461658037973;3.9181763749862450474;-0.060595685135806001687;-0.12172609538556408704];
  IW1_1 = [-0.22136801509961170908 0.20050943844020144291 0.20091471024436632908 0.29046359749142508511 0.20131056413625092461 0.18650551773852289128 -0.42509101176131935951 -0.47382386317273106346 0.15455273783450140801 -0.20627246006625624131 0.47218420239870523281 -0.17350373856547465778 -0.16300302380094935928 -0.51162266536381695303 -0.49817310895324851128 0.021763397561281337678 -0.4574758023664393658 0.36377498731391827436 0.084097404237450698594 0.30098705263086916739;0.488692945525220801 -3.5603631614979738274 0.76787676240003832628 4.5051503183540226871 0.69405834444311520492 1.1980942702732666483 2.2357427749105078796 1.7501259724636391013 -4.6310489850819145374 1.7499080026585043246 3.7579662521097683126 -1.3318588255876300241 -0.82380648550739921099 0.71547841544867440966 0.9197211051574548657 -1.9475773951781434334 -0.71236021764659640532 1.3397475236661953168 -1.0034781035669666505 0.38642559230749679022;-3.0542640134403233887 5.9909298061092011167 2.4667706511208824516 -2.0223276376503838137 1.2046854429254745256 -0.76981995646500156472 -1.3665684908691932975 1.5562347182406877089 1.4328200544218134826 -1.1845922625972111764 -0.081007104315211464485 0.96906400415013571603 3.46747138018951917 0.21208977034532891515 0.53032333609478432557 -0.31959531847312472319 1.0914292260928375455 -0.65246135860155574093 0.92824853110736627038 -0.16616385820189638611;-1.3083420074553335777 -4.220132237306930989 -0.40625376195982093597 -0.17443342686751189818 -1.3229993784107825228 -1.3842798993871914259 2.786396393156199025 -0.50115637168555826619 0.89287194927925450649 4.5125590866993618988 0.64773583257754296927 3.4542350308637006506 1.7114724024100975708 1.5860759231441221573 -0.425535807085454576 -1.2777862369178953905 2.8809857615007934406 -2.8024912238220611727 0.97782635549289986798 -0.36453923401823673522;-2.1180648445100427324 -2.9678246288884460036 -4.6644913744736644645 -1.1735292949660767192 -1.3136765141946322721 -1.1814116514419998882 2.3809119031421359658 -5.0915023334504194708 2.9848839054238780477 6.0438740525608771748 -2.0818993252646391667 3.0435419568541131063 -4.5864180855296714512 0.66799279678977363517 -1.9431784838105741198 2.5529298810683451748 -1.5909813311877840647 -1.0650138030790399224 -0.96001116611212700125 1.7477014349091384293];
  
  % Layer 2
  b2 = -3.6491349841625466688;
  LW2_1 = [-1.6182559103576379389 9.5613235884171867696 -11.194684226592192644 7.2583223037405826261 11.077593224562027885];
  
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