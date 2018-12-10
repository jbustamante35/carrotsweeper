function varargout = ginput_gui(varargin)
%GINPUT_GUI Provide ginput-like feature on GUI and more.
%   PIXELDATA_CELL = GINPUT_GUI(I) inputs image data I and lets user
%   select pixel values from the plot of I. The pixel values are outputted
%   as class cell.
%
%   Class Support
%   -------------
%   The input image can be uint8. Though not tested, it is assumed to
%   work well with these formats too - uint16, int16, double, single and
%   logical.
%
%   The only output is a cell.
%
%   Feedback / Bugs / How-this-helped / How-this-sucked /
%   How-this-could_be_improved / Anything about it are MOST welcome.
%
%   See also IMPIXEL_FIGTOOLS.
%
%   Platform: MATLAB R2011B
%
%   Divakar Roy   2012

warning off;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ginput_gui_OpeningFcn, ...
    'gui_OutputFcn',  @ginput_gui_OutputFcn, ...
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

% --- Executes just before ginput_gui is made visible.
function ginput_gui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.input = varargin;
% Update handles structure
guidata(hObject, handles);

global data_struct

data_struct.gui_ready = false;
%% Backstory: data_struct.gui_ready was added at the last moment, as in
%% cases where the mouse was moved around during the GUI loading, it was
%% seen executing the function - detect_mouse_hover, even before the plots
%% were plotted and the error msg reported was that handles.uipanel1 wasn't
%% found, suggesting that it was executing the section1 of function
%% ginput_gui_OutputFcn.

return;

% --- Outputs from this function are returned to the command line.
function varargout = ginput_gui_OutputFcn(hObject, eventdata, handles)

global data_struct

%% 1. IF AXES IS PARENTED BY UIPANEL, XLABEL AND YLABELS WON'T BE SHOWN.
%% FOR THIS REASON, WE ARE REQUIRED TO DEASSCOCIATE THIS CHILD-PARENT
%% RELATIONSHIP. THIS IS DONE BY REMOVING THE ALREADY EXISTING UIPANEL AND
%% PUTTING A NEW ONE AT THE SAME POSITION. THIS STEP COULD BE AVOIDED, IF
%% DURING THE GUI CREATION, THE AXES IS CREATED BEFORE THE UIPANEL CREATION.

if isequal(get(handles.uipanel1,'Children'),handles.axes1)
    set(handles.axes1,'Parent',gcf);
    uipanel_orgpos = get(handles.uipanel1,'Position');
    set(handles.axes1,'Position',uipanel_orgpos);
    delete(handles.uipanel1);
    handles.uipanel1 = uipanel('Title','','FontSize',12,'BackgroundColor','white','Position',uipanel_orgpos);
end

%% 2. INITIAL GUI SETUP
%% Get the inputs and store with the declared global variable
data_struct.input = handles.input;

%% Control the GUI components with data
data_struct.handles = handles;

%% Status of Buttons on toolbar settings
data_struct.pressed_tool1=0;
data_struct.pressed_tool2=0;
data_struct.pressed_tool3=0;
data_struct.pressed_tool4=0;
data_struct.pressed_tool5=0;

%% Save original pointer type to be used in mouse hovering
data_struct.original_pointer_type = get(gcf,'Pointer');

%% 3. INITIAL SETUP FOR PLOT
%% Store the clicked data
data_struct.plot1_data_clicked = [];
%% Store handles of the printed plus signs (to be used for deletion of
%% those signs after the finish of point-clicking)
data_struct.handles_textarray = [];

%% 3.1.1. Setup figrue data
data_struct.plot1_figure_data = cell2mat(data_struct.input(1));

%% 3.1.2. Get the x and y ranges, so that proper ginput-type mouse pointer could be shown on the image section of the figure window.
data_struct.plot1_xrange = [1 size(data_struct.plot1_figure_data,2)];
data_struct.plot1_yrange = [1 size(data_struct.plot1_figure_data,1)];

%% 3.1.3. Plot Data
set(data_struct.handles.uipanel1,'Visible','off');
axes(data_struct.handles.axes1);
imagesc(data_struct.plot1_xrange,data_struct.plot1_yrange,data_struct.plot1_figure_data);
axis off;
axis equal;

data_struct.gui_ready = true;

%% 4. Initial setup is done. After this point on, the mouse pointer will change based on the axes and
%% their uipanels and also point selection will get activated
%% 4.1. Start Mouse Hovering Detection
set(gcf,'WindowButtonMotionFcn',@detect_mouse_hover);

%% 4.2. To be used till the mouse has been clicked for the number of points to be clicked for
uiwait(gcf);

%% 4.3. Outout the X-Y data
varargout{1} = {data_struct.plot1_data_clicked};

return;

function detect_mouse_hover(handles,object, eventdata)

global data_struct

if ~data_struct.gui_ready
    return;
end

if data_struct.pressed_tool1==1 || data_struct.pressed_tool2==1 || data_struct.pressed_tool3==1 || data_struct.pressed_tool4==1 || data_struct.pressed_tool5==1
    return;
end

estimate_new_position(data_struct.handles.uipanel1,'uipanel1');

cp = get(gcf,'CurrentPoint');
cpa = get(gca,'CurrentPoint');

edup = data_struct.uipanel_estimated_position.uipanel1;
cond1_plot1 = cp(1,1)>edup(1) && cp(1,1)<edup(1)+edup(3) && cp(1,2)>edup(2) && cp(1,2)<edup(2)+edup(4);
cond2_plot1 = cpa(1,1)>data_struct.plot1_xrange(1) && cpa(1,1)<data_struct.plot1_xrange(2) && cpa(1,2)>data_struct.plot1_yrange(1) && cpa(1,2)<data_struct.plot1_yrange(2);

if cond1_plot1
    axes(data_struct.handles.axes1);
    if cond2_plot1
        set(gcf,'Pointer','cross');
        set(gcf,'WindowButtonDownFcn',@detect_mouse_press_print_plus);
    end
else
    set(gcf,'Pointer',data_struct.original_pointer_type);
    % If doing nothing, you need to tell it to do nothing, or else the other condition for 'WindowButtonDownFcn' will be exceuted from the
    % previous mouse position and that is not desirable
    set(gcf,'WindowButtonDownFcn',@detect_mouse_press_do_nothing);
end

return;

function detect_mouse_press_print_plus(handles,object, eventdata)

global data_struct

if ~data_struct.gui_ready
    return;
end

cp = get(gca,'CurrentPoint');
%% Print a plus on the clicked point, to be used as an indicator
h_text = text(cp(1,1),cp(1,2),'+','FontSize',10,'FontWeight','bold');
data_struct.handles_textarray = [data_struct.handles_textarray h_text];

%% Store the image data
img_data = data_struct.plot1_figure_data(round(cp(1,2)),round(cp(1,1)),:);
img_data = img_data(:)';
data_struct.plot1_data_clicked = [data_struct.plot1_data_clicked ; img_data];

%% Procedure to follow when user double clicks on the plot, indicating that the GUI's work is over
if size(data_struct.plot1_data_clicked,1)>1
    if data_struct.plot1_data_clicked(end-1,:) == data_struct.plot1_data_clicked(end,:)
        data_struct.plot1_data_clicked(end,:)=[];
        data_struct.gui_ready = false;
        set(gcf,'Pointer',data_struct.original_pointer_type);
        delete(data_struct.handles_textarray);
        uiresume(gcf);  % If user wants to close the figure window after point clicking is over, please use delete(gcf)
    end
end

return;

function detect_mouse_press_do_nothing(handles,object, eventdata)
return;

function estimate_new_position(h1,handle_string)
%% This function is used to get the new position of the uipanel, to be used for GUIs that are Resizable.

global data_struct

gui_position_org = get(gcf,'Position');
uipanel_position_org = get(h1,'Position');

est_pos(1) = gui_position_org(3)*uipanel_position_org(1);
est_pos(2) = gui_position_org(4)*uipanel_position_org(2);
est_pos(3) = gui_position_org(3)*uipanel_position_org(3);
est_pos(4) = gui_position_org(4)*uipanel_position_org(4);

data_struct.uipanel_estimated_position.(handle_string) = est_pos;
return;

%% Figure tools that are not enabled with ginput and therefore cause problem with 'Mouse Press Down' and 'Mouse Pointer Movement'
%% are to be tracked for pressed on or off and based on these statuses, the mouse hovering has to be decided.
function uitoggletool1_OffCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool1=0;
return;

function uitoggletool1_OnCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool1=1;
return;

function uitoggletool2_OffCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool2=0;
return;

function uitoggletool2_OnCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool2=1;
return;

function uitoggletool3_OffCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool3=0;
return;

function uitoggletool3_OnCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool3=1;
%set(gcf,'Pointer','hand');
return;

function uitoggletool4_OffCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool4=0;
return;

function uitoggletool4_OnCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool4=1;
return;

function uitoggletool5_OffCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool5=0;
return;

function uitoggletool5_OnCallback(hObject, eventdata, handles)
global data_struct
data_struct.pressed_tool5=1;
return;
