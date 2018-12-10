function [Y,Xf,Af] = tpb329d92d_dd49_4ee9_b143_30ab597811f9(X,~,~)
%TPB329D92D_DD49_4EE9_B143_30AB597811F9 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 04-Sep-2017 02:01:03.
% 
% [Y] = tpb329d92d_dd49_4ee9_b143_30ab597811f9(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-839.604414481266;-7316.6928526506;-873.097427232646;-1646.05375716082;-736.048566587917;-346.884767278008;-251.161393954614;-168.317352166511;-342.308065266741;-543.135073587429;-739.08292924996;-243.822193480643;-167.382252836748;-227.78543702633;-166.270203945366;-156.623146275686;-201.347350487774;-110.880020914278;-106.365322333964;-451.498000640036];
  x1_step1_gain = [0.000435258359298379;0.000254266760874674;0.00034735300303884;0.000777513986910233;0.00172505617054034;0.00206247464333953;0.00302249934878053;0.00323098252238339;0.00220614635584283;0.00242812082773066;0.0023189737723531;0.00407495867726386;0.00341201206501257;0.00551378415971643;0.00648970220781273;0.00620026189110775;0.00633602455501979;0.00774864116843269;0.00243421857993303;0.00313801600663361];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.0520671706201030204;-0.37175340945692297856;-0.83849673296281623802;1.4267871304398507259;0.926988841775980843];
  IW1_1 = [-0.17213006003908906671 -0.3868668062045451661 1.0510746158833548236 -1.1811161793215920746 0.21024703968358257611 -0.17999435761218254837 -0.091596926493085956134 -0.050197223556189739246 -0.54103224412882311434 -0.25677842023768621127 0.88349687499085605147 0.46923228549215578731 -0.14938255851054063617 -0.21481408951536309915 0.15504719391410728968 0.10667157942638250123 0.64729610404056336925 0.2035150686333151715 -0.4787650591930866506 -0.13514369122769379072;-1.8443060666867931996 -1.0009005061915818846 0.72218384420320691497 -1.501488567658597173 0.37895784761678069641 0.55977571347665489032 0.21066170965559785633 0.18140653718448171428 -0.60299967622805694667 -0.13264481030299113407 0.53183873796954894697 0.84063099898819726441 -0.66823612368232498948 0.030948138340968288362 -0.54525333143171383909 0.34727314138796405718 0.52442461261095441039 0.88597147778365670234 -0.76529722785153730946 0.13660639465607868348;-5.963611852783432532 -1.9836994853944061035 -0.35003299204052079618 -2.0217991362335308025 1.3672054250763339134 0.2569252538075618153 -1.0725145327532676909 1.6765525465702397767 0.70301461058494996159 -1.007579424013230085 -0.93920864624787070429 1.8352694001454128436 0.72070165269031738475 -1.0420092740204092152 -0.61410319403310753561 0.17411466095157937528 0.82058974565327658102 1.2056494801709980447 -0.033229176385432997176 0.33913168048706926472;3.4722181836345265005 0.98838612156588256941 0.37497131634876157413 1.0426341726209105154 -0.26752560852386342294 -0.11528535657361002609 0.15737538914087692565 0.027356319195298687413 -0.13416587016393438647 0.57804538579327824355 0.39943276443353809002 -0.83424866010026810237 -0.10821660260852179747 0.31673815925578596708 0.59466075246959371192 -0.2254365388610843346 0.10125846006736756444 -0.52900295965382126973 -0.9023296526171811216 0.70583242954852043471;-4.1335482417895645924 -2.9331405123407168212 0.84819022395980514517 -2.7170098783320688618 0.91645959751176464536 0.56849133900340109715 -1.0555082825054931117 1.1984328978238567576 -0.1368421071431755176 -0.52616612798617223845 -0.23905978952750789324 2.2083214409057392302 0.3156725200120110153 -1.0166279366386838934 -0.43185417868781250661 0.30278220044156262114 0.67683046158060000685 1.5114569455015511856 -0.28548046242115582283 -0.095513551848732697325];
  
  % Layer 2
  b2 = -3.9335365145825784694;
  LW2_1 = [-3.2225904762997950037 -4.2341295063978225954 -7.2194022425386128461 2.8910181984333793714 -8.5449905248857263018];
  
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
