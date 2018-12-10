function [Y,Xf,Af] = fftNN5(X,~,~)
%FFTNN5 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 22-Mar-2017 12:56:02.
% 
% [Y] = fftNN5(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timsteps
%   Each X{1,ts} = 13xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 3xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

  % ===== NEURAL NETWORK CONSTANTS =====
  
  % Input 1
  x1_step1_xoffset = [-0.747006382377769;-0.697742834656546;-0.490797845384199;-49.0053039056943;-64.3612323214913;-54.4179107076702;-0.000809726559050569;-0.0013602606530604;-0.00146328627525471;-0.00149716143554126;-0.00168341702480117;-0.00172499056143884;-0.0017998233406095];
  x1_step1_gain = [1.20478154368001;1.48846702877638;2.11403974925536;0.0124437578985053;0.0144281183648748;0.0162992707192427;1036.50533353783;712.60385144272;620.397298639487;686.558993024048;573.223015277061;572.320767918306;560.820488988983];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.4905873430849347905;-1.7425488170475660343;-1.2428791956928995077;-1.0780341649050932951;1.8184540713397674416;-0.78639293368306961707;-0.58948519730415827045;0.82811873434796323679;-0.332882635320568665;-0.14956563717170759387;-0.054138820548383717601;-0.28545030917469205178;-0.32610445650679287244;-0.70691275449897517102;-0.69860984205194986973;-1.0451650912262209214;-1.535976379021294802;-1.334789126199830811;-1.4328315162866511745;3.2833849054079022878];
  IW1_1 = [-0.91419403754344918944 0.71965946132331859442 2.0480291682439752954 -1.3434507944692732906 2.597521055869225215 1.6040383721878881396 -1.2862448185509363086 -0.31713813064576967804 2.7381072711733178693 2.8241404098121152266 -0.33783227190714937294 -0.51167754974094070697 0.39340660693731954733;2.1034459959497882053 -0.60520348833047965531 -0.99067928620849121657 -0.49042142346591210966 1.0304254957664711334 0.34215885358135267325 0.1079346995384564134 -0.39543022795324239693 -0.033306293144624757319 0.76757379100724898446 0.25688077556153710024 1.1170610982018194246 -0.69400147835188163992;0.91766420024503547737 -1.0823903465092876708 0.03400336715940811827 0.26212034662354988246 -0.050180820582575690059 -0.20723287858785513493 0.5391100638206722806 0.52164950977562352108 1.2040680041349978158 -0.25914399451859537749 -1.5293715245222683219 1.0332690062474523263 -1.2455289057946785647;1.2009715745096138928 -0.19393001942622337119 0.18681767548933406187 -0.38543834951264149824 0.65157801406066484518 0.43336986921193465516 1.1558017638763804502 -0.46657362225001097755 -0.35266274871415004277 -2.2197099738080163789 -0.51614102004135031798 -0.46866948103471556264 0.8448111417898908071;-0.40226289101391593528 0.088654589598546135187 -1.9000433295825138735 -0.087543553309262875373 -0.39997015925702911199 -0.42039199208554900444 0.44551701248854019832 0.38713436583344440489 -0.049407927651636332489 -0.152610548106082905 -1.9264876193575637409 0.90020575679539704961 -0.035521235200955299238;0.30499512154414348819 0.0094909198210360017772 1.5559060368767745342 -0.80717250288702313554 0.72894733828037061585 -0.52066732093882239774 -0.10702511391617298975 -0.02627350190535945873 -1.5340733784951137153 0.35486234721798298741 -1.6905526608285681434 0.32656044849528775842 -0.51180844486316434683;1.0385105157084915106 -0.85154045521819732301 -1.1698579229374717059 0.24609972303115201786 -0.42460957368657492594 0.71309642381560478519 1.1776076798308106053 -0.0528860315306858772 0.44595730271472372452 0.057555726661101890007 -0.72260214648541798432 0.39737312698248872245 1.5459646672999662886;-0.53630446939901765813 -0.2357346822586678381 -0.3097164988026308019 0.00074379786293756789978 0.16438481673305724295 -0.15208501714073402455 0.079867356665207334121 1.1372882421105916517 0.52578082851474550186 0.41499610566104488951 0.26544714004314728939 -0.53909253456172312458 -0.0071703660641876578591;1.1080538389993255244 -0.26814010224951900785 -0.19784317871454157323 -0.57459634218997446276 0.21040398704350779524 0.32152775723751814141 -0.087860110514918626778 0.8521714552406368659 0.036084892451771356792 0.26003684685100525442 -0.29866144913822534335 1.7684882181548984281 0.44865635768333533839;1.1164715959222362152 0.42867886431844742789 -0.18399900698139390331 -0.68652161464202787755 1.6734970769782591127 -0.34407529182455615091 -0.19630573205816889693 0.22427980084683216133 -0.43735572273271694499 -0.77847874839991149098 0.72922501132291028814 0.48912057087714372328 0.10840914382522216786;-0.87337305864417613677 -0.28225966035921334907 0.51037898980309814867 0.16461277139853194851 -0.72409346397410090379 -0.036374487646681913366 -1.2342557273836096687 -1.2923155308455855916 -0.10676649194179049784 -0.62772373307195372671 -0.46117282310421636726 0.09403546124892463165 1.1508356170483495795;1.2356376584694275245 -2.7352642501854340651 1.7034159077261883386 -0.018215511186927324516 0.78435213749964505237 0.79423498350805987744 1.833218245210219921 0.17527696021139813065 -0.48763960038004167119 -0.56645354510711076745 0.1013492413807472875 0.032823142609632628697 0.5376118965611537126;0.49821839463801637082 -1.1141513181040496594 -0.88132848925920559413 -0.015923835392376364806 -0.035926108683668002419 0.38666116067974887205 0.85632726185729235446 0.80785233121845412274 0.42673883612354163164 2.1325351198252215923 0.23895378665160238452 -0.88932961271929977798 1.0526683540020369545;-0.8770592720448581181 -0.61693581748857395652 -0.13229364150441089154 0.30901307250672604976 0.31161974199338243752 0.057461633320727109608 0.42137454095687709721 -0.0073346115110098329681 0.41979847540752690538 0.037773938874458116988 -0.12952807759880161087 0.84216139015430391268 0.73841892177196022384;-0.27097542780670547113 -0.52508321219695230919 0.58147510381283795677 0.12370981326701874226 0.024511583180539399385 0.058121043408562111932 0.58466896089595754571 -0.10838570230798347249 0.7450278029369195032 -0.77538872246401369104 -0.10305158370314955463 -0.27079514790409425773 -0.34877398958328115208;-0.0005001777300510357703 0.98065972815112412153 -0.92001526964725410807 0.48940436037827850235 0.23144335368778232165 0.29132458203617395132 -0.25893360876248261437 0.61429557098547060701 -0.12204685730546382949 -0.3883941797573720045 -0.037264796421938362059 -2.6000287422839791773 1.7736081045552081203;1.9695326986919836632 -0.29099965775693908965 0.55607878657439968606 1.4047113519321408059 0.044857137199191238208 -0.21283410017102491474 2.2958054455610996669 0.13698614013940324519 -1.221487056214558109 -2.7867251827338166414 0.9856388327495061441 1.0726502110868134654 -0.35774644480622119769;-0.30757992965166619515 -0.87285199837222327002 -0.074679301711381373963 -1.3296882551523439009 1.1819244296858111287 0.86639301360858944623 0.57670673247122827831 -0.69720345949303175193 -2.3489870095070308764 -3.747557002668088888 0.22198501300375811374 -1.8175370227412768198 -0.85043699054554389161;-0.24358052870339874652 -0.52036337415086009184 -1.1713038399909800891 -0.82730966379472747896 0.98779004816545901146 0.93073807257308460716 0.69541528885604864119 1.4076931505874625827 -0.38200381800106919528 2.7357252643110556711 0.1457327544861279911 -1.7496762038228204084 3.2923126102843593976;3.3231709484655120868 1.7817708079683720612 -1.8608152827170028054 1.3161177782555495419 0.034076111024574341768 0.91868540914776186668 1.6729691591389912464 0.40929087899213084656 3.2223907801568287823 0.56614608527009191441 -0.6839668524149653539 1.3546016821917834516 0.11262166632079122219];
  
  % Layer 2
  b2 = [1.0596651884913592667;-0.30031577865444170072;-0.46809746814141023341];
  LW2_1 = [0.80511982386273095447 -0.0082785711755419960634 -0.64891590270358023851 1.0005760857140553188 1.2191067160309736472 0.88087523908833953445 0.35682072016721205765 2.900299616811049308 -0.54250880345908092384 0.13413270025028178889 0.96281374773879657702 -0.95087616300170285655 -1.2198467414645766116 -0.91587181704471543053 -1.2528502140849080426 -0.12382212118300464643 1.4969959572013131766 0.81862931564000973594 0.21421047836183002011 -0.56009476506948119301;-0.16551309806761504984 0.24999679281847164702 -0.11128590423416775956 0.41073242276831456454 0.11422132192807349371 -0.049011621881741580797 -0.42795484507307873789 -0.077098630264993944783 0.93646720923995496744 0.73667707623038503062 -0.11212679282062608388 -0.12518899481435594168 -0.13931960522230205091 0.96031445359215950042 -0.80256248668555218995 0.7305819856995050543 -0.32199591410146849446 -0.16023766324876734246 0.6151645705214483506 0.75833743998812497278;0.21353944925970810842 -0.9438759939992217296 1.1727263503131690747 -0.42512921756779437077 -1.1272397325126561896 -0.99213454917759191787 1.0531584325489606879 -0.24189731568895286862 0.60888720294098297448 2.383354269172522244 0.13743287669543713392 0.20142653203959975938 1.6798747138456950889 0.94350016612496478974 1.1562914663935044413 -0.69285202052980909837 0.07054124088616996191 0.75527945599872703397 0.053703895231310562475 0.16520179744772853225];
  
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
