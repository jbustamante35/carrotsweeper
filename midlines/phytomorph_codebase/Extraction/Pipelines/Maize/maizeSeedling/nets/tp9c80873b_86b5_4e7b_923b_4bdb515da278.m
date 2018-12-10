function [Y,Xf,Af] = tp9c80873b_86b5_4e7b_923b_4bdb515da278(X,~,~)
%TP9C80873B_86B5_4E7B_923B_4BDB515DA278 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 31-Aug-2017 17:43:44.
% 
% [Y] = tp9c80873b_86b5_4e7b_923b_4bdb515da278(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.91534726369;-869.204807020144;-780.471896742374;-460.205555216351;-234.31060487846;-540.291178237133;-572.142715533789;-192.819908210792;-263.80049063085;-210.820779921778;-258.183122033525;-286.590370604325;-190.32078243269;-193.663740648593;-194.470619647339;-225.820325675023;-213.719449945752;-258.858847398802;-135.174215342509];
  x1_step1_gain = [0.00034254044244578;0.00027748074658588;0.000806013241167019;0.00156677380853905;0.00168393926281428;0.00269046527889069;0.00198084992045411;0.00236361816544284;0.00245784730269614;0.00195235680232595;0.00282106579415897;0.00424658079482961;0.00317261769260757;0.00555458687461897;0.00429724736478088;0.00484270636687457;0.00333931608099047;0.00598913897545116;0.00384783002604038;0.00642488839732228];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-4.6159871585020484375;-5.2132970707741348448;3.8231353563342729274;-0.55547043384545446987;0.5428403126663414735;-0.88544240763830139862;-1.9201287870904677479;-0.26167154121104940767;-0.77351919009397529514;7.6145962018957726869];
  IW1_1 = [5.1431717695161864867 -2.1509312920198024877 1.6556132845553364064 -1.6808156256110828863 2.9556137421512795882 -2.8529692142821683198 -0.9247128807564422015 4.1535971947272027904 -4.0116300343186086508 1.309744630627523776 0.99339583791206753638 3.5546053043051211162 -1.6290652689650431117 0.59888145075795373717 -4.3316019282103566113 0.80389923673624763545 1.6451421942574826485 2.8260934129119714697 2.1892969670639730317 1.4250777041274846901;-6.9223347087098545671 8.6497678937575059877 -0.041881029407108145979 -1.8984813745979081112 0.64479982717452410945 -1.5353759305270582303 -1.3846232629218306887 -3.5973439084026939128 1.6580194632676881206 -3.8594352225086376862 1.675766608982514061 -0.076728199058303969693 -1.6437208302924568581 1.5317246387234202398 1.8847946071971128923 -2.4894466223069668942 -1.0472584201265575654 1.3703549059362367668 1.276583701999395748 -0.29112568600805083241;2.4038980405083050229 5.5012861296936206301 -2.5367567023953818683 1.949881742530860107 7.568621799987688803 1.4972312610537530375 -0.1184169252988497012 -0.61722440581545034988 -1.1981975588527182808 3.7335474443011467116 4.0626223442875009084 -2.4913659260595042433 0.96986422702908503535 -1.1246668725177488835 2.1103252268855405838 -3.5318469947877129123 1.7993872654609019968 3.7927263811919620373 -0.12413527146615445373 -0.30386394788978721149;-1.9022262736616921686 6.3043552040719719187 0.89278619586290053345 -0.85998731032559116993 0.53032789957611548015 2.1600016791698330998 -1.1636189692970437815 -1.1216996275352881796 -0.90590610681181038544 -0.65453979840614440189 0.94179105132680185353 -1.678624215554838317 1.2091364219640345912 -0.60376697389030953289 0.3126048903087571551 -2.5616575562515766507 -3.7196791786270630986 -1.0045337667086757794 3.5537255118123112219 2.2775333164292779209;-2.0918485790649103961 -1.4925045758966062337 1.8667128543601216872 -2.7537429697428010833 1.7570937625306017882 3.0389975095451040588 0.77691570237783547093 -0.56362255895814761253 -0.16846359048572642769 2.6871382743052527609 1.6258758263753052553 -0.20155476109372894022 -2.4404191297285140472 -1.1732625139013510562 -0.42695696940784511053 -0.6808447918536177923 -2.6620454957884533975 1.9461177910905476107 -2.3480734688429154566 0.42395393613707677005;-2.0760754689077898227 -6.302917345785724379 0.42299411714661577788 -2.2538286106405682752 -4.5898184265963948292 2.0129323949568767027 1.2745216134588392976 3.2823557048819993476 -1.5960141432079331114 3.550691549923464585 -6.8074339761461279963 0.6815588341727845112 4.422951859620418702 1.8787138597077914781 -2.0365350673681215099 -3.1896092103329523049 -0.76039225480900407561 -0.2006348236989502487 -0.63876821188605314994 -1.0242460887088038479;2.8011978450580965472 -3.5598225702383436975 -0.38397268965453701384 -0.94179898281747453925 -1.695535415380089006 1.2375120996705042398 0.12547579687462057096 -0.4137736141900741349 1.5117394379847188812 0.76898813796232357021 -0.48247172696878959375 1.4714479815631200665 0.41844696068265657685 -0.65641518398683640356 0.25717518161727881409 -0.54945404006114373363 -0.13429118233267722804 -0.36445414558510436098 -0.29181643286256681025 0.64866095777508325781;1.3808612275179463769 0.6022484702861296002 -2.5856885705360190819 -0.50834760593919747773 -2.5544421515896229558 -2.7415483807783420644 -2.2118021552919637784 -0.2004468629197509455 3.5995059642837476943 1.2494893642843534209 -2.2042778412961014922 1.9082958867299768535 5.2453282543494639256 -0.11164130984242125189 4.4580174692145568827 -1.34850789410124694 4.1615740060692418822 2.0802941722732097851 -0.88134768185491507708 2.5496060253489312153;6.4077227818528159631 0.4386909594631954179 -0.43761422379050202158 -4.4027210771973255277 -1.3471548387667124924 3.7168587972033519939 -1.1996491291397817225 6.8414091008070530009 -3.2434943414553636387 1.3593719521471319123 -3.7491528394490170051 2.0334486106806446415 3.3429337322361609708 -0.12471049676124865146 -1.8257221896902902092 1.2057345326663042595 -1.6667061396828943831 2.3003272642341849874 1.9192880009916546413 -0.029237034853129097167;0.49588673644931163809 -7.6360804699127138306 -0.034056974859430495794 -1.6842687160415532599 0.63369666065652641151 0.66080504408824392026 1.7665654315129639684 -2.7503621222875600516 1.9746550268852189891 2.0784254342879266986 -3.5714863899313646911 0.16744559180305312873 0.85374165907032595335 -0.67313307233959662756 2.7831074426632622654 0.2350679639464709203 0.14723617660126486761 -1.2407891026362942721 -1.3687988630107776888 -1.8635777029516666925];
  
  % Layer 2
  b2 = -1.0546460461033888567;
  LW2_1 = [3.8672492018808215519 2.6014544372614554213 3.2953004656326299049 7.5036661173131173896 -4.9366932220100450124 -2.0440978987772826514 3.6148709340557623193 -3.2505534506078719303 -2.5443727127516133635 -5.0187507876931620743];
  
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