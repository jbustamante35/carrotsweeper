function [Y,Xf,Af] = tpdc5b52cb_098b_469e_b0d7_821c5c7fb483(X,~,~)
%TPDC5B52CB_098B_469E_B0D7_821C5C7FB483 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 30-Aug-2017 10:26:53.
% 
% [Y] = tpdc5b52cb_098b_469e_b0d7_821c5c7fb483(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-5466.61789079662;-862.937878923586;-1028.71862143181;-451.710971542891;-692.538815519995;-441.043793232208;-783.667639824724;-311.950474474658;-656.5558280693;-416.411813428012;-256.518565057556;-268.74339823522;-384.146123020641;-242.574763141859;-313.061223281976;-248.555985162473;-296.419908591291;-275.821003918872;-174.982344511206;-176.062817024311];
  x1_step1_gain = [0.000192953812201724;0.000635170113819406;0.0010737117105794;0.000992967389065707;0.00140247416929523;0.00179426803546528;0.00152094087103249;0.0026247520165371;0.00199415494766504;0.00251075910029925;0.00454433874041204;0.00378662000738737;0.00302237737948991;0.00448879604138586;0.0036608073346032;0.00363057230533704;0.00381924645546824;0.00427537582274444;0.00592120551781387;0.0060242213864429];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [1.2657279515559829441;-1.0986970716641297763;1.5879057366220110659;0.74756036561051053546;0.43969627548066009259;-1.2891717764349914788;-0.91856195155078856907;2.7289013812319957708;0.59764206266196384654;1.0005308997193878007];
  IW1_1 = [-2.7356477617464149965 0.53240036303951132624 0.38724670323183046916 0.084273266158895368538 -0.11846429325288014567 0.49362153653329771119 -0.076097990534285966935 0.0051345471002936307614 -0.21508066067043715708 0.57693298820538074789 0.18153510694080704502 -0.9008338648674744098 0.86301822723675902793 -0.3350530080490116136 -0.44749372140648652207 -0.3379023969883141465 0.015163013988889495823 0.00038071173076839083527 -0.14134676996013839645 -0.055794647026810112456;0.14698115574624107116 0.88006419579273076348 -0.33101085631415660959 0.24347093485820361924 0.74667818452038836607 -0.69512357222696552483 0.52541417165727821725 0.75029575465277065227 -0.75275783837472476367 -0.55446687408458594781 0.26461107471800138224 0.084876099766401966185 0.80314683894162530731 -0.47817368053341030576 0.24847360186194500375 -0.13230039507588631231 -0.29444911139195029115 -0.16552061954372285224 -0.37835004243832742743 -0.41583349599495128546;-0.26461101179109691728 1.3965555116381691558 -0.75326704109779596408 -0.23192850619825536307 -0.1977286728601669874 0.16451714947613993445 0.67596944648130619093 -0.4990791576518341377 -0.96854400340687052395 -1.0584337971276851498 -1.1533600984984708759 -0.019097580308402473115 -0.419074487512833882 -0.11966734719744982762 -0.47664824236621172737 0.47705082490290823261 -0.5380686384909526554 0.79823016132604462936 0.54339725340277933441 -0.40401899107027722247;0.2335538231559840805 -0.49452679992689135835 0.16469771113776193738 -0.54098166862035212965 1.2186581078700822189 0.45092962062565561387 -1.0841538631881659782 0.54814998620961230547 0.56417413060382803902 -0.0098133836443375799946 1.6234076917286601738 -0.4774125298510254467 0.109665948970829083 1.1275014007248151948 -0.21317346197682615627 0.23625546511493067392 -0.45741005285209346276 -0.18598801226720054758 0.15677684173596201367 -0.11060065531316987109;2.6402980831401010597 -0.05354426733385429793 -2.6860783612150780364 -2.2181555389141820456 3.2668207670047286761 1.5312604450773625864 -1.1155333798244491117 -0.56268594192379761676 -3.4706679438067462584 -1.3279512747159172292 -2.7691496489106470769 1.4005800392458935821 -2.4139103657971374695 -1.3487904700769977051 0.67491929297339336991 -0.54859744749499306238 0.28205352205786443642 -1.5444679570599340046 -0.10574104246667465068 -0.20508271479888290356;-0.097819211968324326634 -0.90645882693682611198 0.37796496393501138877 -0.11397208652214330105 -1.2055470993333960994 -0.5137386317579191175 -0.74790132455228619524 -0.13082820990563195918 0.52062366571878393451 1.0507683417031772333 -0.090049052355115344337 0.15648155737824767653 0.3637764391886066373 0.23015392500843168122 -0.67624008861127471448 -0.58303926143871476206 -0.070252395734753003698 -0.22805564474507414308 -0.65022417432191048636 0.42838293875292676738;0.66846867863695369483 -3.0297499677509867055 -2.1467866913612430046 -2.972254172526376248 1.5787505429877939367 1.0481204513179127602 0.81055282018521646048 -0.50183312889137854373 -2.1319371664769786179 -0.77897430437966974637 -1.6383694244035476029 0.67442210459555584912 -0.34797353776076750398 -0.10707339682774071055 -0.67868097380756486903 -0.42399530379845767891 -1.0246254414416269274 -0.066948737367263283504 0.33422839959554462297 -0.45802358611794069532;-0.071584074577934531103 0.65375529409740018583 0.047093223908938122713 0.12388426949510712716 0.23014463856943653819 -0.16602150319679115364 0.25713159532821211473 0.1809539121279767071 -0.23933364940721585379 0.405136542915191189 0.6278101561962956545 0.05561071027886470286 -0.2722606697971213241 0.480456828609022224 -0.38488542691018390318 -0.29085504359490677029 -0.15978908239894104582 0.95199787737233754203 -0.0096191285658987629725 0.11241694818210727635;2.1912212167513542838 -2.0672880242018862873 -2.3577740747508069141 -1.768242324095503637 2.560372728223911043 -0.61873605774778173139 -0.77268583387948885921 -0.4736195843163359398 -1.1221101777611655059 -1.1643967001181414478 0.25936675927932656416 0.11891179317897279666 -0.037980665905576777019 -0.11152895167448864666 0.83692965975186983663 0.56146380558252773074 -0.58178642118354040313 -1.1856716631706067933 0.30151600723209209809 0.45526414677506238382;4.0298260790817312937 2.2570411426983723047 0.27569165805234463473 0.035180110006358594754 -0.57860575034476824374 0.050983191123321988658 0.83279443799678620763 0.43797646901023684629 1.1020809226084244958 -0.083239497065758186189 0.18622613096365930874 0.21412183728047001763 0.14750448934321228056 0.32022865784915016674 -0.16203069363243338685 0.67454179724607321589 0.19623932989203388133 -1.369689002727953131 0.32086599075309429896 1.0949143122634654723];
  
  % Layer 2
  b2 = -2.1185287113370301348;
  LW2_1 = [2.1500446975833829022 2.1905810269254479827 -2.1903559114314870016 -2.405928113995410289 -6.2179685522964369682 1.4125302617098176583 4.4261174985585514108 -2.2260396323221951498 3.0958751298530313711 4.1554512173299400501];
  
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
