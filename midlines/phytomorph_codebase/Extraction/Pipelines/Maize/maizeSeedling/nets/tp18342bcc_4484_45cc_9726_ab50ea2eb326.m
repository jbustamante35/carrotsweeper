function [Y,Xf,Af] = tp18342bcc_4484_45cc_9726_ab50ea2eb326(X,~,~)
%TP18342BCC_4484_45CC_9726_AB50EA2EB326 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 02:54:54.
% 
% [Y] = tp18342bcc_4484_45cc_9726_ab50ea2eb326(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-2001.76879548539;-2528.40271623793;-1416.92748927618;-923.044887924526;-742.627689023618;-590.342170311172;-636.877630731672;-474.328050496928;-408.756773605251;-607.267828431685;-597.158358369107;-252.737018064957;-357.74330150126;-455.326587536772;-262.500890231167;-229.281637810806;-168.878309979175;-299.148398757005;-374.938438931833;-382.085749121029];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296744;0.000975455655635303;0.00136867839470428;0.00177765334299387;0.00147664952758211;0.0019652698447852;0.00211681540578918;0.00191395299404994;0.00213803411115011;0.00313745816598814;0.00239499394158855;0.00262802245354383;0.00340299377226972;0.00307580218448474;0.00420342200465314;0.00353937526659748;0.00305696135696254;0.00320422591704717];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.5486381400167337219;0.74039382209977910865;4.3679769617405339588;-1.0926650137524736017;-1.8514083897675042145];
  IW1_1 = [1.1100845599181419221 -1.0990731202091814644 -1.2170775398304685133 0.45431884326348864755 0.96688843357771181442 1.0550910557051451644 -0.12247362850455970062 -0.94461155808558094904 -1.006695313289337701 0.56859117912468692957 0.80134924260104911653 -0.5333755589859227042 -2.4279148797411247429 1.282377207332493052 0.19772825660756448984 0.052691648881193924037 -0.70262288207413670094 0.23597338037983894843 0.22628698781223363645 -0.73328580408410226621;0.085935866636968630261 1.9064653292265005824 -4.4369550504878878527 -0.39727813112657062167 -4.2241188714576614061 1.2016791707661589328 -3.3541949025836736986 -4.3837498189930341752 1.2991572541150775422 5.9209478633935752612 -3.1187721014508920092 -1.2224377736443177156 -6.025075058759695068 -0.29956343418780778665 -1.9050543920298783629 -0.50396589591332885405 -0.49633391084944772365 -1.006519913779931219 0.57703899529551017089 0.10421845358690957906;0.85310768463685004992 -6.9364872762479388157 1.8703367257105263377 -3.5288469767240977504 1.8325422332995437991 -0.23833642751879047239 2.9739566537637398547 1.0773899588497370949 -1.3617070470478791488 -4.1356246321609537731 1.1466018553823957848 2.2582179542733982203 1.0010992129684090912 0.45121220953701540735 1.6705350360442539071 1.479937328157298948 1.5801146020442509155 2.422321869621435031 -2.6199774658529983995 -0.75282318894744437721;-1.6190179672788715948 4.2125096316606098767 -2.366349221761660182 0.84355447945381722974 -0.28938674298206329416 -1.4106081751053771711 0.26277913487182086305 -2.354335943657602126 -1.5900958260345219308 2.9571405023513261945 4.5512202915837827177 -1.4526465095237690583 -3.0084866600537383263 3.0433059868114051838 1.9356708701199001332 3.1928788939081478304 0.91447169565912200806 2.0171667665075929676 -1.1559124145087789071 -1.7281395341471832428;-2.1461337474969721661 -1.0389321693142627989 -3.4321866437779657844 -3.4860706665826373118 -2.4697977343530852501 0.55798283351170196998 1.6246792598881261149 -2.7073215660961631457 1.8761778648228160105 -0.22738602881052305538 -0.7645173650097201179 -1.8550555207757337683 -0.12371606775230528652 -1.8273727484536619325 0.40935115187319565111 -0.75511249952385006701 -0.81434271880505293417 -1.5667789965826426357 2.2076270323869833234 -0.10212258635831912568];
  
  % Layer 2
  b2 = -4.6152340919131455976;
  LW2_1 = [5.541446083133446443 9.8840259526185523242 -11.237693373190381863 6.0555368096399950772 -8.3855150292351545716];
  
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
