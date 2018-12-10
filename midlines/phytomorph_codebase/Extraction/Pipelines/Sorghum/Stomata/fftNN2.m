function [Y,Xf,Af] = fftNN2(X,~,~)
%FFTNN2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 21-Mar-2017 18:16:59.
% 
% [Y] = fftNN2(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timsteps
%   Each X{1,ts} = 17xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

  % ===== NEURAL NETWORK CONSTANTS =====
  
  % Input 1
  x1_step1_xoffset = [-0.947070134767602;-0.670169622772847;-0.494539180366921;-0.331867816992961;-0.243894504037108;-45.0586813362772;-57.7491605918741;-58.8634582863445;-54.6160586304585;-49.1413743041316;-0.00210691056798095;-0.0034626342052105;-0.0022140798547634;-0.00392027568001063;-0.00350177609670473;-0.00290825072364119;-0.00325775690404614];
  x1_step1_gain = [1.18296228058876;1.5125230178388;2.27467231745216;3.12905141166457;4.43634593172235;0.0128679981859847;0.0153734858379096;0.0165559328189731;0.0188438341787172;0.0199277744006371;308.587456541411;376.736396435924;377.865181452176;291.419439339357;337.615436227973;377.401563596564;313.780260936842];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.2544085329575520582;1.7145353554530062823;-1.1476669233045018181;0.35567303959899299048;-1.1968077927071993916;-1.226494493218779791;0.39438871679320097741;-0.068240591707197525984;-0.076397336286858449905;0.081594069276905753507;0.019931125543642281422;-0.31269816835111491349;0.89994930153092200786;-1.6274060947502131391;0.84004035765247819789;-1.212641870397124233;1.0951228221582063327;0.86010988206080096763;-1.6082786518537761644;2.0431264022500745803];
  IW1_1 = [3.8212552882040786706 0.10550946261619979782 1.7800991576703153463 -1.3243357098061643295 3.1315143191899275088 0.25869357427495115465 -1.5850943220637045972 0.16318105732802085095 1.962236725528532677 1.2632258224209580444 0.31551760506555032171 0.73798997908709185189 0.45030283964772055416 0.41831116069877172237 -0.65833457653057936287 0.74393117581921652803 -1.690221588942813602;-0.35433026682294849197 -1.0060812692190661544 -2.0604871951068282065 3.4930593316210964971 0.088154722109383978901 -0.11957944824198041278 -1.4420752572974204764 -1.1267621094764466161 0.38426475426010769931 1.0053915726754469873 -0.71104259791003754376 -0.18365122593979368859 -0.17584100164499155694 -1.1743572874603951384 -0.6550340886866512502 0.60512420400522282105 -0.40461813318453410115;0.31923320615605388362 0.070099451658177075464 -0.31216031180547226809 -0.2151374014886944519 0.58553334929407985232 0.5192715284447673163 -0.72915177095325567347 0.53045887958145598873 -0.11621288709503541092 -0.36918875566522485476 -0.70368218241557256576 0.25185103424136612293 0.7968211020991452509 0.10862111356505141369 0.72475460280327086426 0.20993601289287336376 0.44608891771343295307;3.0345899114707792776 3.6077365742639218027 -1.8896955803900403748 -0.19212132914284446583 -1.1503153662152512116 -1.1497386005763670092 -0.010014917776195205248 -1.0033424632328695747 -0.044068868268722090342 0.8216932096262002938 -0.8688669492646662329 -0.63917435700286995282 -0.11505262332572133155 -0.90284121429896724287 -0.31863177463698827419 -0.46512194712423388099 0.091208666407786409547;-0.1383418151099063409 -1.1214738832121893264 -0.8797498453234295912 0.47222620877398524541 -0.28290361208983988428 -0.59959211488609143181 -0.25717061205444374794 0.55443418317009229845 -0.080979704522012307955 0.026078495641020067036 -0.29543957522584946007 -0.34302819031660469751 0.36581691763105406157 0.30873640865928492216 -0.8037128360678578165 -0.77566000841227644713 -0.95949363240432172795;0.51261680166398360114 0.29537177289445226691 0.39088735236259514183 -0.85339787662931199108 0.3304652234164467961 0.2494556098034502456 -0.093295625965353790043 -0.085949917041499956927 0.38507968265146996778 -1.2481578686859484684 -0.301327670164892647 0.20416373552722791307 0.63470900669721297405 0.69674166794726444252 0.82714901848086108682 0.13664000822450303452 0.38944512512372747581;-0.81772179593751137539 0.46293565574191192225 -0.4043982497440310464 -0.82034463305434324187 0.45842659605394742162 -0.41504899821306517893 -0.18726989508882091462 0.20241511042564122591 -0.32090419047341822711 0.17353517574345042607 0.82365236622689519397 0.0017519773084287264004 1.6020769847467182245 -0.27389260501989282082 0.32030195808037686955 -0.59010121676346871666 -0.22474592290920902826;1.5952617863608815618 0.63889226992768965996 -0.39077488853638786548 -0.37958601268412694996 3.2422063187499747805 0.30232600226468092641 -0.010106022429181852479 -1.8828235612906885077 -1.3454197282003412806 -0.54180924481770409429 0.399198112081876888 -0.1278893673861892788 -0.50932594604441838726 -0.4013242180779383883 0.52556255568327236993 -0.13692360729781305473 0.14062717113974054772;1.0014278917890919729 1.4127939198548320476 3.0660769847234106145 1.6137948962959371979 0.24842861064156626649 -0.31242981064915997802 0.99507782950826340063 -2.1528073108671397229 0.79606899290877475295 0.56884643172865045901 0.14864442410910239101 1.6090004657493202078 -0.91482760935416662651 0.0552292664550687179 -1.1288448781116431352 0.12103614273408463964 -0.30318805305243384796;1.2562381188818829703 0.30144044667855884878 0.3958799430536118491 0.8865881298303408542 0.28645506320710378834 -0.091103160043278730296 0.20817483031017697259 0.049584890754478291863 -0.77541905645963427052 -1.031110687930924863 -1.5510693637119830068 -0.37190150959273038378 -2.4868538851081347651 0.51789696374187144112 -0.019502928428253660936 0.7227950411649013418 1.7746844881022030282;-3.995174356478215838 -0.32829224114995625428 -2.7049395333170083155 -0.97418343143237462467 -0.79225935024909766025 -0.0090505953024644300586 -0.73207265251760189617 -0.35244497896584126018 0.39780610544722805333 -1.3888085573294690622 0.42289543196559342686 -0.42186517621520119858 -0.35694818379202503111 0.76979336857540370254 0.60606711629242526751 -0.31543284437247165464 0.43999409694209107302;1.128662344773814441 -0.88283560825101470293 -0.39753482158522385204 0.4138172552082719502 0.48941253002528278815 -0.052137579069228993933 -0.093267196794121168257 1.0436149108569261479 0.35850627644005439443 0.57704953826157956431 -0.18315330975138169323 -0.38980088017635866615 0.30558463535569191372 0.29841490741643833529 -0.36531205419229928477 0.053040484925775371883 -0.071928005666761413162;6.2438342769405981159 5.5630672525246493976 0.21982974481683612011 -0.94709170092108807548 3.2481943229719094646 -0.21994058861465817789 -0.18742349383146425512 1.9609273053277542775 -1.084255460544756966 -0.042639554636574684976 -2.2171756074659718649 -1.4288416872673570523 1.7408724445830308092 -1.4991366282350555572 0.59768660174731402801 -0.107623522322397705 -0.44192249208791672954;-1.2058694935356877309 -3.5576788608682647919 1.2791462909748043675 -1.3945945874644012363 2.2220379859496417474 0.59306226578611298539 0.46395628696844259942 1.0443388408175131588 -0.9001041083142049537 0.12813031049505002601 -0.16544490903243982971 -1.9501253394752622317 -0.68707954130499027112 -0.023079414890700720997 -0.34403102678889702082 -0.057499621035421985549 -0.028994140547205175112;-0.296449062179911349 -0.3228886408444671674 -0.14219701050006228193 -0.017703871574209332357 -0.36461791714068592096 0.58942773441600238016 -0.25267678932934800606 -1.3034603483999236406 0.18349330503127164738 -1.1493339326745335427 -0.24302967812089987243 0.30327831255970649771 -0.31373504166151161554 -0.50585153446614095074 0.55565883826652451294 -0.35143406088681355959 0.23867360747970720047;0.56839690932956710157 2.260327163736399303 -0.16516830289989431479 0.93930871327039167085 2.8453633720391744433 0.66370274807907714276 0.81169803737472945304 -1.9799221278720817718 0.98353606210519828323 1.3199200589112394244 0.10310306756881171264 0.97839488486277970125 -1.2051610478743945976 -0.81140996230565332059 -0.71327827103465890612 -0.011337192016931205846 -1.6434972749730172303;-2.4960705200046371388 0.47712444563088562965 -0.96921526753196274306 1.1772603105401693302 -1.7256580232456806279 0.51998080744070740433 -2.6931211255629832202 0.079878682666138789426 0.06762077466846530549 -0.16900736045954509179 0.27955955895361672692 -1.2976105706602734102 0.51509751579380735897 -0.37851021903575160499 -0.1319891233406396136 1.4036704851887666923 -0.86251708216730238465;0.23260720218235705903 1.0368147004445980031 1.3364741467641618744 -1.6273341541182737569 -1.3354905361576225875 -1.2378880130250393687 -0.56240218817329545864 0.27586189508435154938 0.028515345792567807492 -0.92952038407066583936 -0.019683735862647429998 -0.059276308976329239453 0.2361359666027161297 0.13774527374549683234 0.58983439556824579508 -0.051393663892350796574 0.697341558999575728;-0.80110774015474162102 2.7063045144769724004 -4.8199922922275559856 -1.5260450101643172616 -2.5533030504066465127 0.94384723330709785483 0.69875790698233608556 3.3096814181221292372 -1.4948593554794205041 0.36967215608416043038 -1.3854336139815341955 0.63186989799394843725 1.306648847704975358 -0.58318558577228263307 -0.035070304344253538931 0.28806464927722846836 0.37454573595257045771;-0.015780045275876206784 0.94186893520126435675 1.4662922882529405566 -0.3831346097272563811 -0.76510125123593675234 -0.27048631024265640255 -0.059513172026353242106 -0.3048565657602466672 0.45626186137126956144 0.35297029396971152382 0.46835785123738610247 -1.0480238794314318262 -0.22386361573715629625 -0.65256459927595555914 -0.31615360543493159495 -0.40091488299643551763 -0.63664851220983142976];
  
  % Layer 2
  b2 = [0.68740272651803890902;-0.62233756922003735657;0.77484837107276849366];
  LW2_1 = [0.85172334538033023232 0.22552743301458941239 -0.91640262893235335806 0.97695272924526621594 -0.49110814083580067724 -0.94673820963237609316 -0.60585301094122656629 -0.18239985314415183004 0.88148930114675694281 -0.41937034429564107318 -0.58009447924482626657 -0.61173395718305167978 0.4699157205694296513 1.4709214992136003808 0.82034451811852837633 -1.0443951507158453218 -0.83831302541641661286 0.8498253355816287602 -0.12966523420813086398 1.0931178766513867107;-0.61597152900273344223 -0.68273703124080220572 0.042971996029464977485 0.48499813298120686422 0.81742471149393192142 -0.82309627001038265792 0.79985190629423008613 -1.5415202640516425525 -0.074222021769502016797 0.73544417334160183675 -0.84214118978513130553 0.059820946953933640255 -0.085929251109452570168 -0.27969285155559475697 0.51041869802687589441 0.026491323905861671045 0.50893785431031168809 0.85742612333889733911 -0.95640074002842878631 -1.0112358374670793459;-0.10110312916599061006 -0.57121226547444980959 1.9612138099503535482 -0.46410634299347880427 0.29732639209762762311 0.49213391668180944682 0.81520912680928103722 -1.074958770228594096 -0.66899049245504560446 1.1435874652505839677 0.8302969607932613183 0.84775415644349538802 -0.50857321234849250935 -0.41967352983833550173 -1.9199338393450142704 -0.48273610308824677073 -0.47943928727139129009 -0.91952166666262225014 -1.2392573247484619525 0.052155822060188650191];
  
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
