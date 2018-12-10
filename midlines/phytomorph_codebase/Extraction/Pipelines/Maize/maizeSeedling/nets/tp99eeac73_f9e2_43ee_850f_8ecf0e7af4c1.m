function [Y,Xf,Af] = tp99eeac73_f9e2_43ee_850f_8ecf0e7af4c1(X,~,~)
%TP99EEAC73_F9E2_43EE_850F_8ECF0E7AF4C1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 05-Sep-2017 03:02:23.
% 
% [Y] = tp99eeac73_f9e2_43ee_850f_8ecf0e7af4c1(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1660.06224852336;-1969.79757456081;-1377.69545894784;-1183.76616946853;-711.642931835562;-672.522601966801;-447.37811863832;-450.686587111725;-682.76547757345;-606.219316949402;-622.812481590523;-410.03933433275;-413.434624609482;-429.961560001108;-321.305329962883;-359.072122656304;-277.037336897486;-301.078315473933;-400.794724948374;-285.384637074665];
  x1_step1_gain = [0.000244911356138662;0.000507837665104201;0.000859080031338269;0.00102992699402406;0.00116476174985821;0.00153608271588201;0.0018649088281987;0.00229564725234531;0.00140550898923534;0.0018787225220634;0.00209870879134675;0.00237421637966933;0.00249683918257886;0.0027179566948212;0.00287765136194792;0.0032353616095219;0.00314965135524174;0.0039955055636223;0.00295376821702655;0.00375051214875982];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-4.1805642096684625386;0.43926933746012264415;0.030456828849403313336;0.2195502721898323073;-2.6060498672570902023];
  IW1_1 = [-1.6860480943997917525 -3.6892599107058940966 4.8563438999642087879 2.3618782569622407408 0.67733775292820919134 -0.43937763349501129939 0.60328483273720046043 -3.1155819883836404216 -5.84156837420007502 2.8454857559095421138 -2.7608112405480147444 -1.1006872581798750055 -0.19913683835071399608 0.86515250199381554808 -1.1742591391835572878 2.303272500452727467 0.18155475863932724057 2.1797350097846552863 -1.50034313659839591 0.97237589747821107089;-1.3561622135183424298 -0.090509007791077622507 -0.06679460659981008297 0.4219237881048596428 0.31505173678600623122 -0.43614681325175880255 -0.63758188381299663039 -0.50417554632645600332 -0.78181908705275715565 0.19768698257232217297 0.62407954169562718238 0.6326892595249976603 0.23459189154448598091 0.57981177219153068325 0.29518614866870868108 0.27043911146703336623 -1.1126342141986353962 0.84033705958228110955 -0.10307797446170771649 -0.26616910196462256755;1.2435595642875834077 -5.3836658180861549639 0.74740059378260792489 3.5748458495262989132 -0.5851333519814375661 2.5415870381546166179 -0.23268125578930040631 2.3161357480402466891 2.7678762000162038248 1.5407411761788978311 -0.77770816785892304956 -1.310626345622679878 -2.9211159580083254639 -5.0770816564573699736 1.0953959897048464178 -0.53338618626705702752 0.70647130949964564994 0.028017654875715981844 0.85338735870827664431 -0.64878664285256337774;-0.6629763622023516767 -3.3338287094991101434 0.12618389423079409695 1.5882971913744654557 -0.40146984483852549142 2.2176599114687367553 0.50245649036237993723 1.3486389536945713186 3.1272991629380344492 -0.567660886261054487 -1.1163347844203417303 -0.39951889365220155659 -1.3700874249521119985 -1.0404408020401518797 0.34473288747341990224 2.4993363094026808113 -2.1099070206984333886 0.39365965204446862202 -2.204735910265478438 -0.0072200228865967634381;-2.9488039719579788311 -2.9340681320088890516 2.5782916630291423665 3.6218960120600689478 -0.70568881348325807057 0.93943336815308464693 -0.81831098057430617931 -1.1885290022810379718 2.7654546194407383375 -0.70520459836866877801 -5.3482288750625492924 -2.0101397266644482897 -0.59901154926768773734 0.28475339774453461983 -0.043184316321077637069 1.0640121931419732615 -0.36554822008837839498 0.024932475761499827555 -1.1365705145550157873 0.70204978077176727336];
  
  % Layer 2
  b2 = -1.2138066028384706385;
  LW2_1 = [-5.7709865405499272129 2.4171049836496782959 -6.0644191553193822486 -7.2371931347394742318 7.9425985169805377595];
  
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