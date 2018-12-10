function [Y,Xf,Af] = tpd68d6fd9_8870_47f3_8562_cfb9f2e0dc60(X,~,~)
%TPD68D6FD9_8870_47F3_8562_CFB9F2E0DC60 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 06-Sep-2017 08:06:26.
% 
% [Y] = tpd68d6fd9_8870_47f3_8562_cfb9f2e0dc60(X,~,~) takes these arguments:
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
  x1_step1_xoffset = [-1874.99930442497;-5454.9153472637;-869.204807020145;-780.471896742375;-460.205555216351;-234.310604878461;-540.291178237134;-572.142715533791;-192.819908210792;-263.800490630851;-210.82077992178;-258.183122033527;-286.590370604324;-169.741999038737;-193.663740648595;-218.521589353925;-373.104709120354;-120.218367860416;-258.858847398806;-135.174215342513];
  x1_step1_gain = [0.00034254044244578;0.00027748074658588;0.000806013241167018;0.00156677380853904;0.00168393926281428;0.00269046527889069;0.0019808499204541;0.00236361816544283;0.00245784730269614;0.00195235680232595;0.00282106579415897;0.00424658079482958;0.00317261769260759;0.00555458687461901;0.00429724736478084;0.0048427063668746;0.00333931608099053;0.00598913897545098;0.0038478300260404;0.00642488839732218];
  x1_step1_ymin = -1;
  
  % Layer 1
  b1 = [-1.0067968144024115329;0.51404014634490313718;-3.2221327420166669953;-1.0474057282378719247;6.3609361443707799211];
  IW1_1 = [0.605494443525598669 1.3087242024585883371 1.2161685185480379801 -1.7266717022539814774 -0.49894238651947198093 2.9385118381318937608 -0.71692923951028009188 -2.7043474250728665176 -1.0931165690699311543 -3.7553617877278560044 0.93864420150258032205 -1.8195175799602250066 -0.0031491199792332714297 2.4391232353044438597 3.1586461775165877164 -1.7604542229268513864 -0.41818565810855090614 -1.061323442054761923 1.5725918138179892924 1.9526158342585415717;0.55115036599790023786 -4.414691624488223276 2.6804952700623161377 2.1214622155677536064 0.46832620916871836592 -1.9922428232944930926 1.3966732289677010748 -2.7012752582242964827 1.3944818727174022044 -0.35734500837343902147 -0.38652259111042963635 1.411300486434832413 -0.49276384007814572064 -0.17244925467390748164 1.9349972625154696626 2.5084017978512305191 1.4181028498875007937 -0.31385064437268272997 -2.1000680776892171053 -0.60568626642097889157;-0.75582920730228009276 -6.5797116310169139197 0.33358888541030368158 -1.2621277603452842886 -2.658178776833414414 -0.019791476808124416897 0.66330375067147873125 3.971334179233050321 -2.4673753038139984994 -1.2876470717808807809 -1.7672280008230940584 0.51869014790298151318 1.0700286662216551559 1.8071414204327562736 -0.44142436335934970293 -2.0940922042432679184 1.64851195988999355 1.3255958207507223534 0.17640303227805251285 -0.95491825086150106117;2.9396456733515550219 5.4242447428274047638 -0.37213694556977167105 -1.6915681158564357389 1.2152983209081291704 -0.5213721811113978255 -1.9302158350015490296 -0.79207274832927276886 -0.72297297352236766788 0.454201777112960281 1.1025036549666351604 -0.85871276645280580908 0.98421901633982922242 0.14142271233936767882 -1.3601070073801695415 0.078306031743644866072 2.5950061957139176627 0.86633463536041510622 -0.39131243970296053503 -0.9089119103136249489;4.0781922883428975979 2.0603404596872483268 1.832862253384725415 1.0849633386249937228 2.2785171388497431622 1.1351753356361617975 0.072000234394339474031 -3.1271892193971684293 4.0047689340165275595 2.2410495448419105458 1.1675213055730915279 0.93933628490724230264 -1.3761572393333725284 -0.64562526341687487275 -1.1979249228437447439 2.2522284487767207928 0.20768007909228367036 -0.91511951756867826457 -3.3925769013640554306 2.2880310492983539916];
  
  % Layer 2
  b2 = -1.9228793143959521661;
  LW2_1 = [7.8352152896937807824 -7.8164577330173559133 -6.3522695371674897657 8.8697155841606765136 8.7261423385446388323];
  
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