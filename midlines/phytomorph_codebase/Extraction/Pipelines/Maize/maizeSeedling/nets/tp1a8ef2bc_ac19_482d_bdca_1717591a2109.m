function [Y,Xf,Af] = tp1a8ef2bc_ac19_482d_bdca_1717591a2109(X,~,~)
%TP1A8EF2BC_AC19_482D_BDCA_1717591A2109 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 29-Aug-2017 21:39:33.
% 
% [Y] = tp1a8ef2bc_ac19_482d_bdca_1717591a2109(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923585;-1028.71862143181;-451.710971542893;-692.538815520001;-441.043793232209;-783.667639824728;-311.950474474656;-346.375264292298;-416.411813428001;-183.589471391055;-268.743398235222;-384.14612302068;-202.9790070808;-313.06122328196;-248.555985162501;-296.419908591292;-275.821003918768;-162.786711803926;-176.06281702437];
  x1_step1_gain = [0.000192953812201724;0.000635170113819407;0.0010737117105794;0.000992967389065708;0.00140247416929523;0.00179426803546528;0.00152094087103249;0.00262475201653714;0.00199415494766507;0.00251075910029927;0.00454433874041197;0.00378662000738724;0.00302237737948956;0.00448879604138577;0.00366080733460329;0.00363057230533689;0.00381924645546809;0.00427537582274566;0.0059212055178132;0.0060242213864422];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.3612025807311300962;-1.0977021871996586189;-1.8227504069710362966;1.7866374767751105335;0.79770044813118978766;-1.127272668720907367;1.1507558209586690179;1.6990868422629352708;-0.20472218968887623269;1.8161503003879613249];
  IW1_1 = [0.081953385891507246619 -1.2893774418790255876 -0.25692855874350412071 -0.48567244369091616685 -1.2939586531168671968 0.96466656285195462317 -0.11234934058620309194 -0.087245576474226260721 -0.73911289934002966096 0.7194321692631513443 1.218042255039831101 -0.046605176602635939309 -0.80746155874392844254 -0.095259604052906329819 0.084946551274515127994 -0.13746873583819826936 0.44509735689017237981 0.1979322171521609719 0.12873033790642285368 0.11271918512670046908;-3.9476362456271436585 -2.2827014318444112106 -0.014141729804184276914 0.66270651225863752565 0.060713300827888200262 0.047776380224413037001 -0.33752731187755158349 -0.33854158174426668415 0.52320445639149648631 -0.44713401167479943643 -0.41453461921568829363 -0.68476917079715093895 0.50058263057352747616 0.11018544941499448797 0.7712759908164157352 -0.49695245977182456176 -0.76147623849693868436 0.060746754572056471666 0.034217939764769604583 -0.92766826837251414872;-0.019001254183808879983 0.4083731312445564865 -0.13454362854501111801 0.0032443238913169520146 0.22808184164482520062 0.096113723520639793296 -0.498951405855755592 0.734567516374131535 0.29404504387391122711 -0.40674803678270976537 0.40857442219087380941 -0.069662716445996183956 0.35454579431249355626 -0.223153393189587268 0.065601114534697177461 0.41971503424667122584 0.098582744604281249967 -0.31960258181855155657 -0.28697069460618573888 -0.029892274599969303334;0.71568687542408604418 0.089208317728160624327 -0.44357264406624613562 -1.6870019423901583266 0.63961181600784711154 0.12679994218860041189 -0.0076983308583183866075 -0.58377612511365750692 -0.28838901046913900661 -0.21983043773442600965 0.25261837954269994455 0.6487149759928618753 1.112563814463898737 -0.18435364458205300009 -0.60120009940445862107 -0.36263011985507237478 -0.12743228963274239041 0.50048915470969723707 -0.16664746650724204802 0.13484128273570838985;-1.0801847598254412652 1.6915326557084116033 1.8357194625873347604 1.421399090657710973 -3.6624867657394348797 -0.7066569213809031913 0.33173745331019677263 0.51085025552581186936 -2.603465967042153828 0.01361015087502435994 -2.2736624805335492461 -1.543518811136872726 1.5587789288784474895 -0.57795730566196756861 0.93465273084688083038 -0.52678370738280444474 -0.24709270253231171988 0.81385094405247659921 0.29852951948027184326 0.76302296355144383355;0.30821538642497142702 -4.2587198781890105437 -2.7952400013918685673 -2.8708200144753464045 2.0926186095747660865 0.71664489568819644827 -1.2696732437733955212 -1.3189617360353391717 2.6824073876687819507 -0.65469323705464310681 1.1402303973212535304 0.91743986831644142033 -0.59087264271635364654 -0.95675331645918426648 0.9670353825394754077 0.44990692437839341355 -0.98198348778754174226 -0.53941440005621077436 -0.79092185618561838822 0.2425222650163357252;-0.12135591441989695971 0.022265272916695010713 0.29095799708658082183 -1.3872368799625121039 -0.11956006479691969058 0.36898605149389973601 -0.046526173339834747678 -0.74875512335283145315 -0.48364819861016600244 0.79208542635213785577 0.75511220407993706605 -0.46152940057086683234 0.79214150012138262458 0.47469728035528147014 -0.45901718945863173582 -0.35096802053529940757 -0.1432572637198309129 0.40249829945839188383 -1.0488148918816029287 -0.22278086722152992682;0.84666798239582796715 0.95844586438957624086 -2.5315329215690813847 -1.7709250897502819555 3.3788195888198839789 0.61707057404564480407 -1.672938468515281274 -0.57040028735537740001 3.5225593424532197773 -1.6779408986732942122 2.2441410879294605429 1.6706354666658747021 -3.5070534739815522052 -0.057635152682249356226 1.1100112551945646722 -0.012395894799662266275 0.36916313371828468615 -0.27281881533826568553 -0.66498365244802937202 0.085443555450018263131;0.38843037319616252612 -0.31906827861980519012 1.2617636653827315474 -0.2846188027528588349 -0.40241007356737379386 0.41627397466636223422 -1.1556046838148532618 -0.41873403260001612614 -0.56909112447761855957 -0.078350364514307976882 0.40155723676653776177 0.63684574536789262389 -0.43834175589648022298 0.47647322379362977429 -0.32628304741657271348 0.17424373245994545512 -0.56837971864136016542 0.51799890985772734098 0.18058033927660216422 -0.68970013940025276522;0.097833581246673670395 -2.2468908813528316593 0.98125944609341742009 -1.1076659797631531745 -2.5354346172977910712 -0.69519247932663308376 1.235429305288944013 -1.612491445024662795 -1.5395708361242412821 0.1892876808994648774 -1.1382415216410493031 -1.0952477669289355511 1.8499900185032513988 -0.80521510778620875026 0.93954247705354787534 -0.667768301623268834 -0.020363854623932354376 0.17350712393677850121 -0.064810306236227382959 1.0471992035145902822];
  
  % Layer 2
  b2 = -2.8068942390416382615;
  LW2_1 = [-2.4264022513977812068 -3.9295760724374870776 2.1276186502411977841 -1.6475273082588279827 5.016653923725712616 6.8729288024275518865 -3.0207568680055691956 -4.0920463433788079044 -1.8730457370404014483 2.9065436871806986652];
  
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