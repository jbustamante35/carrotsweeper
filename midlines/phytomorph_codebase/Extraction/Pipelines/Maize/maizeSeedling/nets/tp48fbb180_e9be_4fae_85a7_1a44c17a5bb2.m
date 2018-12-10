function [Y,Xf,Af] = tp48fbb180_e9be_4fae_85a7_1a44c17a5bb2(X,~,~)
%TP48FBB180_E9BE_4FAE_85A7_1A44C17A5BB2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 15:35:32.
% 
% [Y] = tp48fbb180_e9be_4fae_85a7_1a44c17a5bb2(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-2001.76879548539;-2528.40271623793;-1416.92748927618;-923.044887924526;-742.627689023617;-590.342170311175;-636.877630731669;-474.328050496929;-408.756773605248;-607.26782843169;-597.15835836911;-252.737018064956;-357.743301501276;-455.326587536784;-262.500890231169;-420.95523700791;-168.878309979168;-299.148398756988;-279.306010526095;-382.085749121045];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296745;0.000975455655635301;0.00136867839470428;0.00177765334299386;0.00147664952758212;0.0019652698447852;0.00211681540578918;0.00191395299404992;0.0021380341111501;0.00313745816598817;0.00239499394158859;0.00262802245354368;0.00340299377226965;0.00307580218448481;0.00420342200465285;0.00353937526659772;0.00305696135696217;0.00320422591704733];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.9127958471583479394;-1.468673391899284697;4.0760346023727063169;-0.53172487602465501322;-1.9041638572369659066];
  IW1_1 = [3.4538679440771766771 -1.2394210553269544572 1.0356358772213127839 0.66162231298001872304 -0.29805498375725408833 -1.6173181371206117873 -0.018518187875004386844 -2.0011280431028977134 0.27432280475459319691 0.74812993583876741255 -2.9714909914503704336 1.2927261057093080154 -6.5038541536510985352 1.7658055792479283586 -4.4830589439448091937 1.1798151416054323271 -0.97493209409603331927 -0.33185209273561694365 0.23105337756883981282 1.9160667646763160921;1.2560732237942389133 2.4716085227956825499 -1.2306047584083257274 1.4755339930911952617 1.7043119485337026031 -1.1554133053504298534 -1.9301975042172965225 1.4194792079894784642 -2.0470609716763643604 3.328876072477028103 1.5116574786716356904 -0.30521719089434773409 -4.1748413252358886893 0.60480442962144187469 1.9053482230586753499 0.69491904297010420599 0.7056728176433053612 -0.84848199603911311772 -0.75033360735390897123 3.6345374816058280132;1.0424454076087346088 -6.1637410140703545736 2.6074292045812397234 -3.3234976533325952808 2.2106180076364596587 -0.10440303881538333386 2.2549055505814830092 2.4131896187604522375 1.1386437870774201997 -2.5416911314922736587 -2.0982960313248155693 0.12138115114976683051 1.6323711973106695972 0.14294467826867290894 -0.3940342609598349588 1.7505374183930157983 1.4852559055281289169 -0.4713432742938123976 3.0571135981536827764 -2.3215982979738791769;-0.36088916071255733309 -3.2488965171814121469 2.1783692424091709583 0.27666984474670941463 3.4000083175682473957 -0.45068975515571213508 2.0598713239911723427 2.9501946677221679849 -1.5023517436464801644 -5.2979933591365284684 -1.707345551088209179 2.8632901432246948836 3.7453238356562152944 -1.6562153640636416707 -0.30603030863680091844 0.59112766412719064402 0.84969417377624645304 -0.29571737158311894866 0.13352675559710588948 5.381969237899965286;2.084745956450575477 -4.1825180560370638361 3.6291711244813260606 -1.3958878063070259135 0.47545635409052683373 -1.0932280755216343682 1.1282947031487953193 1.9779570431449360335 1.3770574690938031104 1.8725300504148716385 7.5241163932932328962 -1.3478041703058818435 1.0856127306772827001 0.17785984734333537549 2.0411520539080005854 2.7655858059220848055 1.936023341891722982 0.43400328954850853069 -0.5764402430847996861 2.1937332311330832724];
  
  % Layer 2
  b2 = -6.3670223384877564499;
  LW2_1 = [5.147507868362440675 8.0486200473965077151 -11.459842217449418555 -9.9470611023035502285 -5.5705294957399829059];
  
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