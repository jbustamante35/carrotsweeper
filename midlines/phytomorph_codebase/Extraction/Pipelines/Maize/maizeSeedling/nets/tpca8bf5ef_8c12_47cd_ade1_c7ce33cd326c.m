function [Y,Xf,Af] = tpca8bf5ef_8c12_47cd_ade1_c7ce33cd326c(X,~,~)
%TPCA8BF5EF_8C12_47CD_ADE1_C7CE33CD326C neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 04:01:50.
% 
% [Y] = tpca8bf5ef_8c12_47cd_ade1_c7ce33cd326c(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-6506.15769134605;-1969.79757456081;-1377.69545894784;-758.119043371398;-1005.44642162532;-672.522601966802;-447.378118638323;-450.686587111709;-682.76547757346;-458.333844314374;-622.812481590513;-432.343869290218;-387.578117402326;-429.96156000111;-321.305329962891;-359.072122656296;-277.037336897485;-301.078315473943;-400.794724948416;-285.384637074672];
  x1_step1_gain = [0.000244911356138662;0.000507837665104201;0.000859080031338268;0.00102992699402406;0.00116476174985821;0.00153608271588201;0.0018649088281987;0.00229564725234543;0.00140550898923532;0.00187872252206335;0.00209870879134682;0.00237421637966937;0.00249683918257879;0.00271795669482119;0.00287765136194786;0.00323536160952195;0.00314965135524193;0.00399550556362211;0.00295376821702629;0.00375051214875979];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [3.9480221721302015681;1.0856218362674723288;3.1380711954600850611;-2.6507457370299074562;-1.1663163712978295017];
  IW1_1 = [-3.7511922528108638808 3.6840412297238867012 -2.5479245116366739943 0.14160691072321882067 -0.19119472667473014216 -1.8116216324297160778 1.8704826196180879272 1.2427210167587747769 -0.69624019804645154164 -1.0903669999620979603 2.026507539571212213 -1.9561891211923141398 -0.43242050955748340124 -1.9690042731832495537 2.3468194628086633813 -2.5706190301629665562 -0.88301868601621513744 0.020165473759856542002 1.0205943456623893795 2.7332616467314001518;0.75914674332580578398 -2.1132437238799086643 -0.63847548926993469287 -0.9421190032544947135 0.58015326146078027847 1.8813114278329508 0.58299865411489260048 0.82898113384755733879 3.9497868462898440001 -0.22810008603611112155 -0.70291558663311570765 0.087157081357765425755 0.3959139348265435232 -1.0115410066207641826 0.051517518158703166919 1.2457659019518487753 1.0478196414002105552 -0.76542990376266695218 -2.9139244869695506424 0.21212944585502258033;0.96451716032185985661 4.7435113011304093078 -3.6441048928831887288 -2.0340012145934025511 4.4772950150063826058 1.0484638015613703121 -0.055462922823701993336 2.4755699708868079689 13.117521581583499568 4.4283773043153162874 0.16355617081178452921 0.66753520743658401049 -0.246736476045023162 0.61897849537765692318 0.22398532279671320988 -4.304216621411452337 4.6485646582964532669 -4.3890758350390042253 -2.8602277290850688374 3.2792639719938776466;3.0865762596941364038 -0.068254290720511118229 -2.3309349366995806285 -2.3382246346051438479 -0.92852321642965929271 1.3427550092347557698 -0.28419562341184750887 2.4029890606046349077 4.9058142518115879227 0.014557965362025179767 1.3074051880601829101 -1.7056038416282688353 4.4069775980909655644 -7.9815779178701005137 1.6314375720506195666 -1.6687326232192221109 3.5300839547319498024 -1.3165223022227279248 0.3426360712762765437 1.066019172310940677;0.78500615990884770046 1.7790602076863635261 -3.8690534841246724262 2.7701156736233607525 -1.0386235973309965175 -0.08303409564730709902 -1.1963424146489542998 -0.90928818994333637793 -0.08617251445602412685 -1.9205164302844406343 0.064671358828397668139 -1.5085044492347183809 0.88535565975292263019 -1.5331412284385097511 -2.0030787634466520331 1.4989100039476872528 -0.61793272152435552158 0.98341720021110812233 0.34818828888986314007 -0.54929943048764373614];
  
  % Layer 2
  b2 = -0.82668015931362903714;
  LW2_1 = [-4.6710515793290783293 -11.892738646967130123 6.2778317018968179397 -9.9694884080876988719 5.7803546629670226054];
  
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
