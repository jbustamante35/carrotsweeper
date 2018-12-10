function [Y,Xf,Af] = tp8b144746_3fe7_4fa0_a6ad_12bfb4f2bb30(X,~,~)
%TP8B144746_3FE7_4FA0_A6AD_12BFB4F2BB30 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 03-Sep-2017 11:06:18.
% 
% [Y] = tp8b144746_3fe7_4fa0_a6ad_12bfb4f2bb30(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-7409.93644730062;-2639.1224350771;-1040.58146532707;-1175.93373113978;-1366.53951203701;-713.224898465728;-638.058802539164;-604.796304556693;-546.419507096588;-435.830146318432;-311.040988170545;-332.822721465588;-407.901071240805;-247.407476323938;-310.614407522934;-284.621358963184;-348.992764796122;-139.23912248536;-228.698126275206;-273.249410813689];
  x1_step1_gain = [0.000226370837385414;0.000563347828528141;0.000948114745185448;0.000858032283205086;0.00073024238319975;0.0015171639926653;0.00146770964020618;0.00151546500149802;0.00206494173697165;0.0020247507057702;0.00326439309746569;0.0022999722951874;0.00248425272472354;0.00412752991038746;0.00322794766513105;0.00343797851616054;0.00328303385537599;0.00723300236088635;0.00445727486740905;0.0040752789910078];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.9350958737342984772;-3.9329015551113646332;1.8464785921234188049;-4.513778982917930449;-1.86482894897483642];
  IW1_1 = [1.7733030290599907097 -0.81253508581245792453 -1.3820637280404750591 -1.0570020810175815296 1.7182677289359244988 -0.93673640679424707844 -1.1355470343659597177 0.42638240355511597413 -0.85892760963647951655 -1.0997378798859673843 -0.84944885468313846566 1.3460204429864424114 -1.4416289297911641931 -0.13418021720850289991 0.20730037376893598045 2.6215197946060393441 3.4712886738598935921 -0.78212996736773732209 -0.45990589914515334069 -1.906464188582152941;2.1601238567519707168 6.3617213093362661169 0.7624518693142928738 -2.6864481316692612189 2.315780208284652808 1.8580427973609465742 -2.1895188807385572716 -1.5668666920677880583 -4.8518373478063949733 0.41372477506651322843 -2.1215047951586050168 1.4954547783695690821 -1.1384553946137505598 -0.38183834055406917773 0.27956029185045566665 0.24932713877016132731 1.8821151157566486933 -1.339586191389925629 -0.8490851023684797827 2.6997695049300269865;1.296038474005294594 -0.97004834062427969243 4.1939220931523086833 0.049901512673050367475 -4.1098338716077380539 -0.27305933425573186701 0.13583405742299364816 0.63618562575238868373 1.5546038012980489817 -3.4537255563445592443 -1.3574090393867945803 3.8969766628694770816 -0.085030099335716813336 0.53312308156123422442 -1.1911393739416833437 0.30449703018057083792 0.63434650245792967294 0.48731270454136943471 -0.12588668467563096653 0.021586413005901260859;4.1111775984894904923 4.4466146018875667423 -1.9654430457207348049 -0.37533525220011182411 4.4936274478004394339 -1.9484009932686603417 3.2958039203959885732 2.9169793079962182958 2.9157056330351851514 -2.9168152941510587794 0.85203642906324328354 0.99623755581902373457 1.319857850632164542 0.68836857254152816754 -0.61656024423536959223 0.10166335978138120455 -1.5731082026031184551 0.041442676158618740878 -1.5200612014984753451 -3.08416104795321111;0.7445138040118225442 -0.98181123389832514192 1.3170284074104803285 -1.6599277922086885795 1.6243633018977905458 -0.36112197973786425154 1.4041286059251079088 -0.22794966468471070109 1.3924248942634604997 -0.85574467418165534038 -1.1519598853250498927 -0.45998166120123545886 0.34623468187260592055 1.1506552875222324506 0.20546191287424686589 0.19801834182141531127 -1.3544508281381129056 -1.3472675163965994738 0.06547945595416344855 1.1234750810446199054];
  
  % Layer 2
  b2 = -4.7845785285153530353;
  LW2_1 = [6.3278909471190392111 -4.296444541366027714 9.6199011110348742193 -7.5936352767762178573 6.2312018303682359388];
  
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