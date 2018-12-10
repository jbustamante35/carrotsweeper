function [Y,Xf,Af] = tpaa87df0c_e724_4229_8d21_d9fdb7eec8f0(X,~,~)
%TPAA87DF0C_E724_4229_8D21_D9FDB7EEC8F0 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 02-Sep-2017 16:58:21.
% 
% [Y] = tpaa87df0c_e724_4229_8d21_d9fdb7eec8f0(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-2001.76879548539;-1846.28602050405;-1096.59977942709;-923.044887924526;-742.627689023618;-590.342170311173;-636.877630731669;-474.328050496927;-408.756773605245;-607.267828431699;-597.15835836911;-384.72168072541;-357.743301501248;-305.70191787155;-262.500890231176;-420.95523700793;-306.924499680699;-265.923075543788;-279.30601052607;-382.085749120988];
  x1_step1_gain = [0.000250772026915415;0.000457175383291266;0.000795694570296745;0.000975455655635302;0.00136867839470428;0.00177765334299387;0.00147664952758212;0.0019652698447852;0.0021168154057892;0.00191395299404991;0.0021380341111501;0.00313745816598808;0.00239499394158863;0.00262802245354391;0.00340299377226956;0.00307580218448474;0.00420342200465305;0.00353937526659758;0.00305696135696244;0.00320422591704722];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [0.6612392405413427543;1.560887888240955812;1.5368059588554641159;1.8180982933640850163;2.6332811261115414148];
  IW1_1 = [-0.2907521163904263628 -2.4867666255897340477 1.1108317342187872256 0.93768020898892423531 -1.1705198068713993376 -3.2243469880259749694 -1.4598161259445316951 -2.2788339854800470086 -1.2509142783926134879 3.8868008849131912008 3.3983658154277640584 0.25453502392211391214 -2.3835456570492632444 -2.5026652111582334825 0.65439447113603177453 -1.6241482656036458909 -3.0446020367998953837 1.734480248518445844 1.4713233123559217308 -1.4326930942328186536;0.19997803888877527245 2.5231062228220295651 -1.7014780495879733291 -1.2494844441600130303 1.4970533845157323327 -1.2707502492431783647 2.2858241990317771375 -0.75883000582675230739 4.242327132246190402 -1.940236584654021712 -2.647638365944162242 0.74292359490340964534 3.4040563305986175635 2.154371768167994361 -0.3684467037436690573 1.2922530281891215775 0.098748436073040332661 2.4183952170045754926 1.8134731999445909967 1.1788871102130824564;-1.2507807199872191806 3.6654869068843916047 -2.5655898806978636451 -1.659846104487215479 1.4152370514232450738 -0.11320047406733992068 1.0641615743697754315 2.3420432623334197508 1.5966826810557617211 -3.8126351410389887064 0.24006861187912037492 -2.5737845875964810993 7.3140315480215543076 1.6305688494687187173 0.80421414263781720866 -1.2405777644418023442 -1.3253309807214823657 -1.3405307812307691862 -0.67224032339550154891 -1.7066191448632028749;-1.0350574074234257793 4.5730865817853683097 -0.54926898282008507213 -0.36549567928942272577 1.2314982416239776963 0.20266322779191320302 1.9633757440167447506 -0.41759393815587214638 -2.0386971504300874791 -3.2411792318976395677 0.18568793322170074056 -1.0324624183191004612 3.9643000856059185111 1.5118854601994238251 1.1917025670474745702 -0.84624421381437986422 0.2224714893051210618 -1.6931542669693651071 1.2441175836351701101 2.8749568708549197993;3.1195826892830527122 8.8902181060938989532 -0.96644672753288707323 -1.6221188009726028589 0.54075743618035154725 0.62485801808421925774 1.46204113518208767 0.54415116367148996801 -1.3549823928244082438 -0.54077002078928670947 0.69566688278748278762 2.0552494224723218963 1.1404523807467024099 -1.2693429214588676235 -0.19974224731080683637 1.5779373516364678665 -0.74396424999603816275 2.7818033450496684189 0.093069686350501834582 -0.25801178676214980845];
  
  % Layer 2
  b2 = -3.2709666095327252222;
  LW2_1 = [7.3830352102678631354 -6.8939891370230563794 -10.839001682717883668 -8.3844211410174285248 -9.4701119240443514258];
  
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