function [Y,Xf,Af] = tpa008c594_f391_411d_928c_e847074b03f1(X,~,~)
%TPA008C594_F391_411D_928C_E847074B03F1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 06:09:32.
% 
% [Y] = tpa008c594_f391_411d_928c_e847074b03f1(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5310.24104946905;-801.674192791829;-1831.68761111049;-751.778745457491;-962.87968250126;-561.234234396124;-567.818719457227;-566.708636588772;-460.813029179949;-650.116847383027;-365.70388877515;-279.141698173319;-445.050021503464;-362.027565171382;-290.108959264272;-350.309427866398;-205.100647540257;-201.92345259531;-223.562425332945;-192.300211130329];
  x1_step1_gain = [0.000187746732766858;0.000517587446146428;0.00076092913976903;0.00107683880161886;0.00113410115771159;0.0015463795383062;0.00191937835803698;0.00192268942981792;0.00211835684506626;0.00196958059572074;0.00268955398662549;0.00257377430529595;0.00284234016821657;0.00279908218986809;0.00409640929201516;0.00345584818973335;0.00422505290367573;0.00549698021899214;0.00471432618261885;0.00476878567270279];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-0.84826727217829955841;1.2726686416048567096;4.6700153299018101194;-4.449937957066078198;-3.1623104163046664716];
  IW1_1 = [1.7241622035899404342 -2.6720218103831214762 1.0687489637468188253 0.6276975368661328325 2.4395653981755307527 0.20862418136795246593 -2.3211017958676816164 1.3807818155606070842 -0.40590739703812300787 -0.36979252418839941852 -0.59472557725360863934 0.042496706253276950171 -0.79897020889180137004 0.75512784239959729327 0.58469407872042078456 -1.0008028399760560134 -0.50590979744919195582 -0.99681298786644068066 -0.8799514520586808608 0.69206841767627658157;7.6490410170878586626 4.8476804192377311864 -0.39690237722041288304 0.59144265587598476852 4.7730945434016733842 3.3678112023849084977 -5.0137996541501443915 -0.81009210506815820807 0.075512022156731292699 2.4151887402309766273 3.2454288260746992911 -0.40262719210408426385 0.8438033420259795081 -2.7262507499085133489 0.97662639046086285877 -1.9941396179552188794 -1.3801342643266438781 -0.46265091286813586624 1.7709890743957801273 -1.9665326140057164395;0.5395084022608764851 6.1148940490150751614 -0.89810113062508134529 -1.1008753570710825276 -0.65081411481960915744 2.7034708228803583374 0.46941746435826320338 -3.0431739856353261509 0.20556008105400050412 0.84102244562008943163 1.6049906947735703699 0.12111156263587348481 0.69422044158432882899 -1.9464940223242042094 -0.71620770440294267978 -0.41978615869170221586 -1.2150020075733785063 0.14503792882053517799 2.3167038532218509062 -2.0054637280990661452;1.0182832619639634597 -3.5005438624286071558 1.0282792686521657899 1.1608821941263849364 1.4498928618298660975 -1.0557309526060649763 -1.2631692345385194809 2.0782503058540666707 0.87468869239418978179 -1.0410650316626477441 -1.1630329154855503226 0.5552203536988281618 -0.86679760691926632798 1.2211684311424029215 1.1088047405098115927 0.36133917408635085611 0.74016868065443908176 -0.13072631921347688255 -1.3428019376861435319 1.4148528216396245583;1.2392412483668888701 -2.4929346782319985465 1.2122365192632360564 0.57292088765553916829 0.71378984359980301999 -0.56051047574940260354 -1.1137803531073862118 0.80783152358562515527 0.27173026462222688027 0.26905612892140651216 -0.60790994684688237637 0.083857502654859278346 0.19408609132882656545 0.56400881312316752947 0.74396522866952785602 -0.36345658130173641442 -0.0066307445130925407248 -0.1746020607055526197 -0.70965696020345370165 0.77328910002658624023];
  
  % Layer 2
  b2 = -1.5810651795737560299;
  LW2_1 = [4.4866355379271682935 -10.836065095719739304 -9.7218352919169372939 6.4628266181492826092 3.991683380921934976];
  
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