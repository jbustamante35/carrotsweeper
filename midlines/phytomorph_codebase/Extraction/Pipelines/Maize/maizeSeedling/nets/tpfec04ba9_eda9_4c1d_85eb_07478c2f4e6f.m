function [Y,Xf,Af] = tpfec04ba9_eda9_4c1d_85eb_07478c2f4e6f(X,~,~)
%TPFEC04BA9_EDA9_4C1D_85EB_07478C2F4E6F neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 21:15:20.
% 
% [Y] = tpfec04ba9_eda9_4c1d_85eb_07478c2f4e6f(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-8184.90793418461;-841.368302509969;-1769.06834334395;-1424.63292253515;-1115.63541033651;-861.953076142274;-996.664725988398;-421.625702789072;-721.553885286555;-557.274532479563;-502.132109545495;-541.944235264821;-543.17700519132;-167.207345824935;-344.528961662821;-294.401776402713;-391.010406370098;-235.62564790017;-174.294824484551;-301.03830724203];
  x1_step1_gain = [0.000202107905561193;0.000554219398814733;0.000742185534904077;0.0010185341647551;0.00103332655312912;0.00147681597144008;0.00118221307527392;0.00265958726683859;0.0012523589192447;0.00165285753178001;0.00232860128237177;0.00188013184297681;0.00198822427969292;0.00479968082327788;0.00252047169225132;0.00363117206158383;0.00272984263620494;0.00383651951608927;0.00502254586493056;0.00330771328744503];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [2.437449340175700474;-1.710965613404849428;-2.8309040974194354945;-0.6146594429918981195;2.5837389660616025822];
  IW1_1 = [-1.2414070291846153182 2.6908815041206306518 4.5634367078172868304 1.3120796625878443642 0.53715623977614701534 -1.7036224669220458861 2.3033759315805633072 -0.32085514551377181114 3.7328434787316044918 -0.47644019357348849075 1.2639801090085562851 0.12145637955930065066 -1.2450971323235699995 -2.9534959626205750105 3.1398621495463179265 -1.7553201346353106693 -3.2927801017795697014 -0.0015761933335365898787 0.91486128169309555158 3.0627522660067394433;4.6852453689927031633 -3.9394617633210748942 -1.1401234524968295592 -0.90914642257889621835 -2.2970113217897893954 7.1444377705392394518 -1.4302833670630246132 -5.581920086698148431 -2.3482876387140945162 0.9399902002394014966 -5.786905995127901825 5.6056516171799728454 0.18955910753582388972 7.2482619273856148112 1.8210060548558504312 0.77221585593126074176 0.039229219183555832928 0.13444074764929531129 0.25571955680211644335 -2.6662948619559085905;4.7148173162909774447 3.1051014801297918133 0.21612003969070933684 -1.1918487843378771096 -0.79600561703939110014 -0.018341608991395932926 1.374600289203117498 1.1665754552868259353 -1.5210311914176402048 -1.8674447566587093483 1.4738514820985233467 0.59499854324434897634 2.2652148454884661888 -1.0640441279454886381 1.0801252878553981862 0.088303251409362196189 -1.99850110979908413 -0.031205317478045030422 1.5443480451480744442 1.3999825138377719469;-2.3641649029303803964 -3.1584330791221657542 -2.3015611297935314461 -0.66312652776984482017 0.53573852752356165841 -1.3481834072581220241 0.65132419331990210054 -0.11677063934030815218 0.50855961087742529436 -1.9885021335626915651 -1.823862217726880397 0.24988715571377462243 -1.8467651804793894232 -3.4325687690011474906 -5.9671479742795927592 -4.4294701155336051457 5.3982207570459044277 -3.390108753759331961 -1.5533960454290645092 2.4751412896620847981;3.3912048876817317833 4.3788372469965040068 -3.4233324465934731329 -1.7535777306914064333 -5.2487780469388276572 -0.0021379691594778853052 1.4146503293405081614 -1.035101242496006213 0.82184730924520965623 -0.79870225722463517037 0.54934595267715435885 2.5017710020137720051 4.1185641521081723937 -0.327040280040294562 2.9022478803674651715 5.6619798296664667347 -4.0134178522522807597 0.86831749081623632591 -0.31380660994965087207 5.0244270001192257169];
  
  % Layer 2
  b2 = -3.1989420953656173552;
  LW2_1 = [-8.4865156763904767701 10.826726058422462273 8.3013768044485658493 9.2216487451669095776 -12.118704557030833158];
  
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