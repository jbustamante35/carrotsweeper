function varargout = modelViewer(varargin)
% MODELVIEWER MATLAB code for modelViewer.fig
%      MODELVIEWER, by itself, creates a new MODELVIEWER or raises the existing
%      singleton*.
%
%      H = MODELVIEWER returns the handle to a new MODELVIEWER or the handle to
%      the existing singleton*.
%
%      MODELVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MODELVIEWER.M with the given input arguments.
%
%      MODELVIEWER('Property','Value',...) creates a new MODELVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before modelViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to modelViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help modelViewer

% Last Modified by GUIDE v2.5 13-Sep-2013 16:11:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @modelViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @modelViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before modelViewer is made visible.
function modelViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to modelViewer (see VARARGIN)

% Choose default command line output for modelViewer
handles.output = hObject;


handles.model_mean = varargin{1};
handles.basis = varargin{2};
handles.sz = varargin{3};

handles.pc1 = varargin{4};
handles.pcE = varargin{5};
handles.var = varargin{6};

handles.XC = varargin{7};
handles.gr = varargin{8};
handles.tipAngle = varargin{9};

handles.frame = varargin{10};

handles.lastVec = [];
handles.xout = [];


handles.modelFigure = figure;
handles.modelAxes = gca;

handles.closestFigure = figure;
handles.closestAxes = gca;
%{
handles.growthRateFigure = figure;
handles.growthRateAxes = gca;

handles.tipAngleFigure = figure;
handles.tipAngleAxes = gca;
%}

handles.histogramFigure = figure;
handles.histogramAxes = gca;

handles.totalFigure = figure;
handles.K_axes = axes;
handles.GR_axes = axes;
handles.T_axes = axes;


plotModel(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes modelViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = modelViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
displayK(handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
displayK(handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
displayK(handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function [] = updatePlotModel(handles)
    mag = 3;
    pos = get(handles.slider1,'Value');
    pos = round(pos*size(handles.pc1,1))+1;
    p1 = get(handles.slider2,'Value');
    p1 = p1 - .5;
    p2 = get(handles.slider3,'Value');
    p2 = p2 - .5;
    p1 = mag*p1*handles.var(1,pos)*handles.pcE(:,1,pos)';
    p2 = mag*p2*handles.var(2,pos)*handles.pcE(:,2,pos)';
    vecT = handles.pc1(pos,:) + p1 + p2;    
    plot3(handles.modelAxes,vecT(1),vecT(2),vecT(3),'k*');
    updateTRTA(vecT,handles);

function [] = plotModel(handles)
    handles.modelFigure = figure(handles.modelFigure);
    plot3(handles.XC(:,1),handles.XC(:,2),handles.XC(:,3),'.','MarkerSize',1);
    hold on;
    mag = 500;
    U = handles.pc1;
    E = handles.pcE;
    for e = 1:10:size(U,1)
        quiver3(U(e,1),U(e,2),U(e,3),E(1,1,e),E(2,1,e),E(3,1,e),mag,'r');
        quiver3(U(e,1),U(e,2),U(e,3),E(1,2,e),E(2,2,e),E(3,2,e),mag,'g');
    end
    
function [] = displayK(handles)
    mag = 2;
    pos = get(handles.slider1,'Value');
    pos = round(pos*size(handles.pc1,1))+1;
    p1 = get(handles.slider2,'Value');
    p1 = p1 - .5;
    p2 = get(handles.slider3,'Value');
    p2 = p2 - .5;
    p1 = mag*p1*handles.var(1,pos)*handles.pcE(:,1,pos)';
    p2 = mag*p2*handles.var(2,pos)*handles.pcE(:,2,pos)';
    vecT = handles.pc1(pos,:) + p1 + p2;
    %plotK(handles,vecT);
    updateTRTA(vecT,handles)
    

function [] = plotK(handles,vecT)    
    K = measureModel(handles.model_mean,handles.basis,vecT,handles.sz);
    mesh(handles.axes1,flipud(-K'));
    axis(handles.axes1,[0 301 0 200]);
    caxis(handles.axes1,[-3*10^-3 13*10^-3]);
    view(handles.axes1,[0 90]);
    updateHistogram(handles,K);
    
    
function [] = updateHistogram(handles,K)
    handles.histogramFigure = figure(handles.histogramFigure);    
    tmp = handles.lastVec;
    xtmp = handles.xout;
    [handles.lastVec handles.xout] = hist(-K(:),200);    
    plot(handles.xout,handles.lastVec);
    hold on
    if ~isempty(tmp)
        plot(xtmp,tmp,'r');
    end
    hold off
    axis([-6*10^-3 20*10^-3 0 3000]);
    
    % Update handles structure
    guidata(handles.output, handles);
    

function [] = updateTRTA(vecT,handles)

    K = measureModel(handles.model_mean,handles.basis,vecT,handles.sz,handles.frame);

    updateHistogram(handles,K);
    
    delta = bsxfun(@minus,handles.XC,vecT);
    delta = sum(delta.*delta,2);
    [J sidx] = sort(delta);
    
    ta = mean(handles.tipAngle(:,sidx(1:10)),2);    
    gr = mean(handles.gr(:,sidx(1:10)),2);
    
    mesh(handles.K_axes,flipud(-K'));
    axis(handles.K_axes,[0 301 0 198]);
    view(handles.K_axes,[0 90]);
    caxis(handles.K_axes,[-3*10^-3 13*10^-3]);
    set(handles.K_axes,'Visible','off');
    
    plot(handles.GR_axes,gr,'k','LineWidth',3);
    axis(handles.GR_axes,[0 301 0 2.5]);
    set(handles.GR_axes,'Color','none');
    set(handles.GR_axes,'YAxisLocation','right');
    
    
    plot(handles.T_axes,ta*180/pi,'m','LineWidth',3);
    axis(handles.T_axes,[0 301 -10 100]);
    set(handles.T_axes,'Color','none');
    
    figure(handles.totalFigure);
    
    
    
