function [Y,Xf,Af] = fftNN(X,~,~)
%FFTNN neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 27-Feb-2017 19:25:39.
% 
% [Y] = fftNN(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timsteps
%   Each X{1,ts} = 12xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

  % ===== NEURAL NETWORK CONSTANTS =====
  
  % Input 1
  x1_step1_xoffset = [-0.999988063784976;-0.999758456590177;-0.999387796587003;-0.998605429935175;-0.998704305628923;-7.71781814732938e-06;-3.25637753246335e-05;-7.46240634303018e-05;-0.000171196095275532;-0.00012397055228582;-0.000130122655381214;-0.000186339763895969];
  x1_step1_gain = [1.00001823825867;1.0001365049056;1.00045518348839;1.00088294899638;1.00309038731579;146623.95883664;41050.9232076262;18185.1235677422;8517.61552517888;9148.82067911263;4860.5675721233;5405.40190460618];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.687033904538931095;-0.35838304700190631857;1.2243127175560926112;1.3004188743138345163;-1.3529015407346627331;-1.1909772064973183081;0.74356133256049283808;-0.196155766691968797;1.0151013039532872373;-0.58043815353487437925;-0.47848392593759853098;0.62241718527700229391;0.52968626976734078937;0.48632845487780879123;-0.87529147513341221121;-0.34697598810030494398;-2.7700763809202042509;-0.69964633447723634241;-0.90251312827258345362;1.2026542069287153858;1.1426963790040185476;1.2485245707769594414;-1.5122950981658465075;-1.8126000591969240183;-1.9720318175458151355];
  IW1_1 = [2.5783075712737746343 0.73055215792612770986 1.8746099807287608385 -0.21054145596492662706 5.6727017758082585175 0.64003479057495826776 0.34397477991250441365 -0.70510691894658661383 0.037619402308918907518 5.107443460055434592 -0.49055936311387643434 -7.7624684101583669715;0.48224481262163265693 1.1592886249227718487 -1.4541141029889208536 0.60716076195699542062 -1.0054480778701924137 -0.87074885786597255866 -0.9603021685827064502 -0.86473023752361122263 -0.58228499315929915792 0.96745180922355722419 -0.7680490103256867318 0.85855493892527701139;-0.40385365088217428653 0.23056567606343106314 -0.65982567386044044255 -1.3910148157074273545 0.3821192732057283048 0.73647413955335083724 -0.0085030501965469729858 -2.2663339668448978514 -1.4624927060340002782 -6.8270041154016132268 0.61909587813541711032 -4.8282707101280166739;-0.09980139578750976681 1.2614233373625807211 -0.98387045143555396098 0.79921502482617612362 0.32119874226630901193 -0.58316286005995277009 -0.25914387808330874297 -0.39244409612681940835 0.48086767334144647768 -0.18072926226628227564 0.038209155666932718975 6.400613155567798529;0.20038336116973040624 -0.15410264244354701657 -0.6940085044023693861 -0.059527977270757606587 1.0125464389672749199 0.49164315099936584108 0.6712435809914354401 -0.23311302148992329997 -0.7105586149911579108 0.74284481220574405391 -0.58810600772168020178 -0.56279608860643637325;0.17676806335250067059 -1.2760828872805152212 -1.2625492482073339406 -0.77050641866181368567 -2.256915894157975444 -0.40216577692589183801 0.53115049239926126479 -0.35112874306094982169 -0.16466063980835976555 -3.6117906094208129808 0.11643336235943535506 4.9019417975792585196;-0.14158827989941902281 0.31686329957410025449 -2.6765613084898030039 -3.3597825001316006244 -0.40823162149684544264 -0.64186979734394289387 -0.71663618882696966672 0.50428409182413513889 0.32060222073454852643 0.0098558975490795203145 -0.70512752531135891587 17.727211303761439609;0.71570029437319859777 0.71050076629312386434 -1.0164466240639822292 0.16548607788592900603 -2.6108014252067084371 -0.78036082144223717361 0.44741298470163853862 0.50374693457144892328 -0.11692130826271777966 1.9265298198496605409 -1.2076583330557220108 1.7762498145634166846;0.06508700651237814272 -0.46309675071330513862 -0.084344926860594984785 0.93232013770624178051 4.8428599872652426939 0.052896027398661553398 0.48080770331478245305 0.60452235122457065319 0.49060174310431653444 -4.8230442970946718617 0.64773352450383070433 -12.805904804498601024;0.079592416525449252718 0.19868647016575044506 -0.72138467955489915528 -0.070211604178992012559 -1.2190143303358955684 0.39657349131889624294 0.21362978408071717218 1.1494807538645399614 0.014651504722334779754 1.2297227288125593958 0.25279551426246260615 4.8010817140949741955;0.018561374460643681461 -0.60874568082123858836 0.33768944213898355322 -0.19074321296197679909 0.091919375127857097518 0.1045481014077825066 0.60124749565196899947 0.81111448496017324317 0.23951531176229978737 -1.127114568348539958 0.99737433956483267927 0.60227809485867633832;0.18848505307459953451 0.16915743455927134686 0.36105491214385332377 -2.6763630613657203305 -2.0219538729684711775 -0.2860075520796961368 -0.54775100774626683808 -0.62995695199732149039 0.66437102637529077409 2.3415631606914906193 -0.28339089439947823257 4.9814740301769306896;0.46756504770507750779 3.5510214484351121378 1.7528642564430887418 1.3212004361665452556 -0.78435651976054687751 0.52516159517739657847 -0.19148428655076690785 -0.47076137180260679838 0.017169956824039087073 0.97745367885775291672 -0.26054006377618399437 8.1106468535973998968;0.29045933920827088359 -0.11863231770067869664 0.010075470290754494698 1.2094309543514794036 -0.74500974620051618746 0.048056432981453331421 -3.4513679325846085 0.26824013406246677294 -1.1261786190870475277 1.7614712129769292126 0.85303067514038688479 0.9072169425798544129;0.1656007137733137291 1.3278678912992549765 0.80086151896518009696 1.7136849500834971938 -3.224649711238247729 0.16395426835840232194 1.5698445447356603744 3.1513658128441983663 0.26156525898499333183 1.6710850939334584986 0.47248265717580867351 9.1808624042200648319;-0.6723961966705302995 0.4996831137974171444 -1.6265037689065948712 -0.61341530550543765798 0.44298807088762898188 0.22760715512164836349 0.25545257142050875609 0.17499301074119080757 -0.015269578179267273566 0.37606995605164794672 0.80957120208942090489 -2.9015342212433141711;0.53362353685729368991 -0.52090935867276100257 2.4574036066180169335 -2.2641166555412080008 -9.8056344608030787668 2.1442761972609876153 5.2161335351524993342 8.5412787614150236237 -0.44711591864251426509 -2.1779132248956774021 -0.72512817762921844 14.542325081691082289;-0.58882331759261141713 -0.0060675125552665919865 0.13374509390356220662 -0.03751800856464806988 0.36197912142926169654 -0.40205259501700113534 0.28990938719682907143 0.6426039992510026444 -0.4172183898069399266 -0.35243592136189888508 0.52093523579717138094 -0.63800147350750735509;-0.894491762483911379 0.59167259268089411606 -0.14935080361190872233 -1.0922185186681134805 0.49339983136623505322 0.41213925558286307327 0.10376170000477062672 0.047359615542284377865 0.67436376500699590775 -2.484474946118099048 0.20366680653731020212 -2.2452168899314779438;0.14266162999952286539 0.80678398487521374438 0.19569707946430583978 -0.18080444692125613448 0.4909719017568786037 -0.5843144812043900016 0.30422104016051826347 -0.13397334890697976117 -0.055136779650182314583 -2.0981126276616559601 0.92584928336071437016 1.2882224744426367113;0.99719655406898388694 -0.60090499256616047852 -0.25141086519676086874 0.074271531085832634411 -0.73698697427301773644 -0.99430509676428380228 -0.027464433003512870451 0.28643296205309537461 0.75339544168196870277 1.4792591400420671555 0.14215712026277363411 0.3457615596540398184;0.93746943949139838193 1.8830140234552654732 -1.4784503564397584174 -0.29023178900448237671 -5.6440765540825124091 0.35719657098930301942 2.1304154167257065566 0.76668356077736798948 -0.34783451657406250224 2.3656498765725393341 0.061607181742859118412 -32.487635034744144491;-1.1666699617851858939 -1.147295931876419095 -0.48136089663599418431 -0.31016821914894188161 -2.6164581755753690828 -0.17455257784704833468 -0.026149341059982587898 -1.1889213324124643378 0.46736292977522891912 0.24827726007771794281 0.6062885995006004336 -1.2657757622285947807;-0.32771300956311588681 0.15135109507070804158 2.2370989435437644843 -3.4407383792566608172 -0.50804362824073312499 -0.21247280039741855973 -0.085246158347971306979 -0.44837981176666136252 -0.3252969894777690274 -1.2072890704754513092 -0.41203417675452252977 -2.9399196853612838964;-1.4050264168950226118 0.17627328150702098464 1.548713053198017997 1.771091685286875661 -1.2679054649767373952 0.65723816230706799679 -4.6596105985757239054 -0.84289260087806316957 0.49000043534794501809 0.14460132952625417979 0.1551717340499283726 0.79067618685981699844];
  
  % Layer 2
  b2 = [0.66204140968060121164;0.51256349037728521356;-0.87251859638253659313];
  LW2_1 = [-0.80297194202446364741 -1.309137612342821555 2.2248743348819139953 -0.48361539669603470504 -0.61924851799213764814 -0.13939255328748245688 -0.61752148949133378597 0.61094796950175889183 0.30396113453067458154 -2.0376781103997996247 -5.2876301154209874866 0.12518374928560699333 0.27816778545725051597 -2.0850441316026455674 -0.15568245182269355165 -0.81571474071196203859 -0.56021707252492713991 -0.53965763610858130672 -2.0400298807366303855 2.3329868805791056374 0.2656297855433157884 -0.55445437221453675036 0.6652623269941116968 -0.060096515637579298041 -0.92641747191107459702;-0.15368001004654352482 0.10299045267291906391 -1.9260769060614224557 -1.0349292978882480565 0.20103794834422997195 0.19546654806288851769 -0.4434356791473944015 -0.5215316133537203358 0.29958716525227657712 1.1333217808000897175 -0.1211174067517348657 -0.25234433740807155067 0.77361192349103091725 0.83510688066745364377 2.34139549306060335 -0.18893999753305851241 -1.0786481077689853425 -0.27813584418593495506 0.49074762794867843319 -2.5419858865654290092 -0.23510943963702696236 0.42622897791677211154 -0.63984501495173307095 0.5376794260837352013 -0.047013687737728161731;0.87542759334241604385 -0.52706562885914298722 -3.7833366121784957059 0.20095892438680282299 1.1021669214498557299 -0.19795014853427284063 0.6413007596183090131 -0.7734131951336179478 -0.44469872856565445529 1.9804722029403085681 8.2046804854569188592 -2.2157876643519940352 1.1430656759397896938 0.77641470838414650135 0.67685865687312474215 0.3812158970009639769 0.4749247627889090495 -0.26997977151138635721 3.0095565492264046981 -0.82752121828183444752 -1.9734785289944118958 1.3695493929737194705 1.0401596470540046546 1.0728803444693328917 0.0036295735819822814647];
  
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
    a2 = softmax_apply(repmat(b2,1,Q) + LW2_1*a1);
    
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

% Competitive Soft Transfer Function
function a = softmax_apply(n)
  nmax = max(n,[],1);
  n = bsxfun(@minus,n,nmax);
  numer = exp(n);
  denom = sum(numer,1); 
  denom(denom == 0) = 1;
  a = bsxfun(@rdivide,numer,denom);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end