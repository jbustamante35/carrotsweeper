function [Y,Xf,Af] = tp520c1702_54aa_4c96_a40d_7efe895005d4(X,~,~)
%TP520C1702_54AA_4C96_A40D_7EFE895005D4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 09:23:32.
% 
% [Y] = tp520c1702_54aa_4c96_a40d_7efe895005d4(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923585;-1028.71862143181;-451.710971542901;-692.538815519995;-441.043793232211;-531.307871775743;-311.950474474654;-656.55582806928;-380.160028036829;-256.518565057549;-268.743398235231;-384.146123020683;-242.574763141968;-313.061223281902;-302.321350914853;-296.41990859129;-191.974128700832;-174.982344511361;-176.062817024285];
  x1_step1_gain = [0.000192953812201724;0.000635170113819407;0.0010737117105794;0.000992967389065703;0.00140247416929523;0.00179426803546528;0.00152094087103248;0.0026247520165371;0.00199415494766512;0.00251075910029925;0.0045443387404119;0.00378662000738718;0.00302237737948954;0.00448879604138452;0.00366080733460377;0.00363057230533666;0.00381924645546919;0.00427537582274555;0.00592120551781178;0.00602422138644417];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [4.014922419322297209;3.5505441281110270424;-1.6228337986594620901;1.5891746116926539489;-1.2299402995507602299];
  IW1_1 = [1.3270040487493228731 -1.7146852553975582278 -3.8125868537796443469 -0.49494427387926487327 6.3924028603407334614 2.3209918968512188719 5.5601354251810954921 -2.1456842704141085854 -7.0651598435664499931 5.4484963979533205247 -6.0150641979569750006 4.8602498170656538434 -8.64670732369682149 -3.13752501154204122 -2.5894088877034140417 0.62125396100707863667 1.3511585771984682136 0.44781960656897068507 -0.54216453869054936909 -2.3079793890892918995;-1.5875405407893272702 7.6421217502887914819 3.1187013609084348609 4.1064965301440548018 -1.2364554057766292683 -0.60798567866137576132 -3.7704980948431434662 1.4429384442748978135 4.4093076919971556649 -1.5737995344721844582 0.95544726243160005286 -1.4142261925267081146 -0.70866506733817014041 -1.4435789211299898671 -1.7745100295072748775 -0.85485024664499986713 -1.0563792975306278699 -1.9741524827465128933 -0.36726327496675009732 -1.5361024094481288138;0.037089570301025928667 -0.10799228950993930432 -0.38393857039240264495 0.83457333807124278735 0.40976337366620002589 0.44682605713143797033 0.06405108796902007906 0.088649431854353158422 -0.28014995738711329087 0.91069353756741633354 0.0029197509334167944146 -0.80278666266005049845 -0.394248350797905045 -0.21873743189841465817 -0.11494882658512604878 0.246568309694131238 -0.60832173974646186476 -0.11783302456775807121 0.17592538729780793361 0.034856051265202289691;1.5596414509059799425 -0.42450229532774358399 -0.84068715151040995792 1.0106024731374434467 0.6534405590994675217 1.2008534205728289557 0.035816300982703691558 0.055449131035503022857 -0.91251204353243231537 0.73855438372782333811 -0.61529337515763371513 0.024239950871983436537 -0.46381570875620192762 -0.29333056166946452903 -0.89450754246503727529 0.51329721075495460259 -0.13660862963884370935 -0.033007665026741357661 -0.38649318667226401081 -0.075951184366238100321;0.92010079567030644832 -0.52927811907846700201 0.40896535045750259663 0.57562454116057404363 -0.82398079852393724121 1.2469851484150225041 -1.6323833639090452596 0.13062489053436859843 -0.77794090975515162167 0.81383588616898983048 0.81884198906233796045 -0.23472120210014033037 2.1364539413553575997 1.5767748259943255018 0.35572484448804730928 -0.17869193344610290186 0.071469768691057145027 0.80486449334741849704 0.86354689752625024557 1.4501596031801564646];
  
  % Layer 2
  b2 = -5.6077215742308688107;
  LW2_1 = [-13.922678938258989945 -14.601844194027023605 5.4265486768952921892 -0.32867910148460538489 4.3641012558969638491];
  
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