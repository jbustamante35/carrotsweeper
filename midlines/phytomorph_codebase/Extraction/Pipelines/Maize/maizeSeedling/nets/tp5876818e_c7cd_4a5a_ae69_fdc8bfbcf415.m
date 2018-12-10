function [Y,Xf,Af] = tp5876818e_c7cd_4a5a_ae69_fdc8bfbcf415(X,~,~)
%TP5876818E_C7CD_4A5A_AE69_FDC8BFBCF415 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 23:17:29.
% 
% [Y] = tp5876818e_c7cd_4a5a_ae69_fdc8bfbcf415(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-2001.76879548539;-2528.40271623793;-1096.59977942709;-923.044887924525;-742.627689023617;-590.342170311174;-636.877630731671;-474.328050496919;-408.756773605263;-607.267828431686;-597.158358369112;-252.737018064957;-357.743301501209;-305.701917871497;-262.500890231166;-229.281637810812;-168.878309979176;-299.14839875699;-279.306010526093;-382.08574912102];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296746;0.000975455655635303;0.00136867839470429;0.00177765334299386;0.00147664952758212;0.00196526984478523;0.00211681540578911;0.00191395299404994;0.00213803411115011;0.00313745816598813;0.00239499394158856;0.00262802245354413;0.00340299377226972;0.00307580218448469;0.00420342200465328;0.00353937526659784;0.00305696135696224;0.00320422591704723];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-3.9227389509795451694;-0.83938435178961823269;2.0505466168812862549;4.2952702042385730863;-0.65837500000496917973];
  IW1_1 = [-4.5184188196706331908 1.606641295440330186 3.5187486678382886574 -0.84769069061375568186 -2.3884670283732729601 3.0037012604708852592 1.3823721312902295111 0.16656821873818999924 4.1860387621946930636 -0.2501945169203343422 2.2899868162961394269 0.54963227215956234062 7.1515821036938866939 -0.36065647314862386841 1.4653506971276293314 0.44350519392984094624 -1.5501787938904931785 1.5734740489358798321 -0.2031093753549757519 -1.286822167381753701;-1.560000743850073146 -4.7335192628511491009 -2.2146152811072408184 -0.9266652830230712734 -1.6924839758707019843 1.8867149644226930416 0.017860617090474407986 0.20273230974151326689 -0.83955487831109687491 -1.0089356134940543619 1.7605692809855382208 1.0907502940653870382 4.7415170811446394339 -0.80403402578090277242 1.1280279916139765994 -2.6626506586175606728 1.6715576264878932999 0.27721737301366200512 -0.49793169312044316399 2.2889949161757670026;-0.77915396709795281538 -3.5196611782885502961 -0.31935002744826157528 -2.7241144400139538639 0.98163193993687491989 -1.050081168902074813 2.2268332089332298906 -0.11765110911075843392 0.81068282446285400145 -2.4095283523287851501 -2.8249454562057314533 0.099738888858756702405 -1.6131651508337145895 1.0892466196271137768 1.9739985921483991227 -0.34537229707169503357 -1.2029918897156610669 1.8221327080778926266 0.98304700563975522254 -0.92784075549391709092;4.3759704156878660442 -1.7969598969327518212 -4.331710898845412494 -2.5078303829246229029 2.7288159311229307491 -1.5808517490357665292 0.11592307427896096639 1.1691254954496661522 -0.47038390115244427259 -2.7056636915326413018 -0.6726815922798111691 0.060088841052051825131 -3.3707908093660936544 -0.23408192206224470722 -0.87821147873651406179 1.8451148796203595825 2.0351231267667229119 -1.6670381764435195926 0.7484328931307375754 1.8198476541836237086;0.26691737283272698456 0.056248620335316780461 -2.6938796766725010556 2.3790223813884301407 2.6357323896707205435 -2.5037713165650727909 0.046184589999293895746 2.0100487014994041779 -2.5196133384981931336 -4.0722603651437765393 3.2073465387120774217 3.3553296819647293603 5.582543447843161033 -0.2656308094814004428 -0.14238031997176620047 1.6667597394058701887 0.16138679077432832587 1.3580027409666972638 0.55526810874487209091 -2.0515982554183098685];
  
  % Layer 2
  b2 = -4.0930882881291710262;
  LW2_1 = [-7.0288611050714697726 -3.3260601593484424576 -5.4609249648942945754 -5.8261600451729229633 -5.7932626003218139488];
  
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