function varargout = iPlant_ver0(varargin)
% IPLANT_VER0 MATLAB code for iPlant_ver0.fig
%      IPLANT_VER0, by itself, creates a new IPLANT_VER0 or raises the existing
%      singleton*.
%
%      H = IPLANT_VER0 returns the handle to a new IPLANT_VER0 or the handle to
%      the existing singleton*.
%
%      IPLANT_VER0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IPLANT_VER0.M with the given input arguments.
%
%      IPLANT_VER0('Property','Value',...) creates a new IPLANT_VER0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before iPlant_ver0_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to iPlant_ver0_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help iPlant_ver0

% Last Modified by GUIDE v2.5 29-Jul-2015 14:03:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @iPlant_ver0_OpeningFcn, ...
                   'gui_OutputFcn',  @iPlant_ver0_OutputFcn, ...
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


% --- Executes just before iPlant_ver0 is made visible.
function iPlant_ver0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to iPlant_ver0 (see VARARGIN)

fprintf(['hello fix']);
%set(groot, 'DefaultFigureRenderer', 'opengl');
set(groot, 'DefaultFigureRenderer', 'painters');

% Choose default command line output for iPlant_ver0
handles.output = hObject;

% selected data for analysis
handles.selectedData = {};
handles.iplantUserName = '';

% mount the private data
[handles.local_irodsMNT,handles.iplantUserName,handles.remote_iRODS,handles.remote_iRODS_return,password] = irodsMount();

% set write permissions on return drive
setWritePermissions(handles.remote_iRODS_return,'phytotest');

% set remote dataset
remote_publicDataset = ['/iplant/home/phytotest/Public/Data/'];
localMntPoint = '~/phiRodsPublic';
mountPublic(localMntPoint,remote_publicDataset);

if nargin > 3
    list = get(handles.methodsList,'String');
    idx = find(strcmp(list,varargin{1}));
    set(handles.methodsList,'Value',idx);
end

%{
if ~isdeployed
    import phytoG.locked.Bpersist.Bfs.implementations.fs.Bfs_irods;
    handles.irodsFS = Bfs_irods();
    handles.irodsFS.accessResource(handles.iplantUserName,password);
    
    import phytoG.locked.Bpersist.Bos.implementations.Ostore_irods;
    import phytoG.locked.BdataObjects.fileSystem.implementations.dataSetCollection;
    store = Ostore_irods();
    store.accessResource();
    res = store.query('DataCollectionName','masterDataCollection');
    % if no master node is found then create and persist
    if res.size() == 0
        masterDataSet = dataSetCollection(store);
        masterDataSet.setProp('DataCollectionName','masterDataCollection');
        masterDataSet.persist();
    else
        masterDataSet = store.get(res.get(0));
        masterDataSet.setOstore(store);
    end
    handles.store = store;
    handles.masterDataSet = masterDataSet;
    
    % populate data collection names
    dataValues = handles.masterDataSet.getDataNames();
    dataValues = jarrayTocell(dataValues);
    set(handles.dataSetListBox,'String',dataValues);
    
    % populate image list collection
    if numel(dataValues) ~= 0
        handles.imageCollection = handles.masterDataSet.getObject(dataValues{1});
        handles.imageCollection.setOstore(handles.store);
        dataValues = handles.imageCollection.getDataNames();
        dataValues = jarrayTocell(dataValues);
        set(handles.imageStackListBox,'String',dataValues);
    end
end
%}




% Update handles structure
guidata(hObject, handles);

    

% UIWAIT makes iPlant_ver0 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = iPlant_ver0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectDataSets.
function selectDataSets_Callback(hObject, eventdata, handles)
% hObject    handle to selectDataSets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selectedData = uipickfiles('FilterSpec',handles.local_irodsMNT);
set(handles.dataListBox,'String',handles.selectedData);
guidata(hObject, handles);


% --- Executes on selection change in dataListBox.
function dataListBox_Callback(hObject, eventdata, handles)
% hObject    handle to dataListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function dataListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in methodsList.
function methodsList_Callback(hObject, eventdata, handles)
% hObject    handle to methodsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function methodsList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methodsList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if no datasets are selected -- select datasets
% get the choosen method - first list and then selection
methodList = get(handles.methodsList,'String');
sidx = get(handles.methodsList,'Value');
% if needed run the files selection list
if isempty(handles.selectedData) && ~strcmp(methodList{sidx},'Train Particle Learner')
    selectDataSets_Callback(hObject,eventdata,handles);
end
% set the selected method
para.jobType = methodList{sidx};
% set the choosen data
para.fileList = get(handles.dataListBox,'String');
% set the outpath
para.basePath = [handles.local_irodsMNT filesep 'phytoMorph'];
% set the iPlant user name
para.iPlantUser = handles.iplantUserName;
% set the localmount location
para.local_irodsMNT = handles.local_irodsMNT;
% generate jobs
%jP = generateJob(para,handles.store);
jP = generateJob(para);
% execute jobs
jP.run();
% transform jobs
%executeJob(jP);


% mount the public data sets
function [] = mountPublic(mntPoint,publicDataSetPath)
    % create irods mount
    cmd = ['icd ' publicDataSetPath '; irodsFs ' mntPoint ' -o max_readahead=0'];
    % make the mount point directory
    %mkdir(mntPoint)
    % mount the drive
    %system(cmd);


% set the write permissions for phytotest
function [] = setWritePermissions(D,U)
    cmd = ['ils -A ' D];
    [o,r] = system(cmd);
    if isempty(strfind(r,'phytotest#iplant:modify object'))
        cmd = ['ichmod -r write ' U ' ' D];
        [o,r] = system(cmd);
        cmd = ['ichmod -r inherit ' D];
        [o,r] = system(cmd);
    end
    
% --- Executes on selection change in dataSetListBox.
function dataSetListBox_Callback(hObject, eventdata, handles)
% hObject    handle to dataSetListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns dataSetListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataSetListBox
handles = updateImageListBox(handles);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function dataSetListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataSetListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addDataCollection_button.
function addDataCollection_button_Callback(hObject, eventdata, handles)
% hObject    handle to addDataCollection_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%import phytoG.locked.BdataObjects.fileSystem.implementations.dataSetCollection;
dataSetName = inputdlg('Please enter data collection Name...','New Data Collection Name',1,{'new_dataSet'});
tmp = dataSetCollection(handles.store);
tmp.persist();
handles.masterDataSet.putObject(dataSetName{1},tmp);
handles.masterDataSet.persist();
dataValues = handles.masterDataSet.getDataNames();
dataValues = jarrayTocell(dataValues);
set(handles.dataSetListBox,'String',dataValues);


% --- Executes on selection change in imageStackListBox.
function imageStackListBox_Callback(hObject, eventdata, handles)
% hObject    handle to imageStackListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imageStackListBox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageStackListBox


% --- Executes during object creation, after setting all properties.
function imageStackListBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageStackListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%import phytoG.locked.BdataObjects.fileSystem.implementations.imageList;
folder_name = uigetdir(handles.local_irodsMNT);
folder_name = [folder_name filesep];
%%%%%%%%%%%%%%%%%%%%%%%%
inPara.language    = 'matlab';                
inPara.filePath    = {folder_name};
inPara.fileList    = {};
inPara.fileExt     = {'tiff','TIF','PNG','png','tif'};
inPara.returnType  = 'set';
inPara.verbose     = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% configure@inPort@indiv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iport = inPort();
iport.sourceType = 'structuredPathScan';
iport.source = inPara;                    
% pull from inPort                    
iport.pull();
iport.fileList.xForm_iRods(handles.iplantUserName);

% get name of new image set
imageSetName = inputdlg('Please enter image collection name...','New Image Collection Name',1,{'new_dataSet'});
%tmp = imageList(handles.store);
%tmp.persist();
tmp = iport.fileList.toJobject(handles.store,handles.irodsFS);
tmp.persist();
% put into current image collection and persist 
handles.imageCollection.putObject(imageSetName{1},tmp);
handles.imageCollection.persist();
% update the image list box - high overhead cost for now
handles = updateImageListBox(handles);



function [] = updateDataSetListBox(handles)
  

function [handles] = updateImageListBox(handles)
    % get the selection from the dataSetListBox
    selection = get(handles.dataSetListBox,'Value');
    selectionChoices = get(handles.dataSetListBox,'String');
    selection = selectionChoices{selection};
    % get the image collection
    handles.imageCollection = handles.masterDataSet.getObject(selection);
    handles.imageCollection.setOstore(handles.store);
    % get the data values
    dataValues = handles.imageCollection.getDataNames();
    dataValues = jarrayTocell(dataValues);
    set(handles.imageStackListBox,'Value',1);
    set(handles.imageStackListBox,'String',dataValues);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over dataSetListBox.
function dataSetListBox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dataSetListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = updateImageListBox(handles);
% Update handles structure
guidata(hObject, handles);


% --- Executes on key press with focus on dataSetListBox and none of its controls.
function dataSetListBox_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to dataSetListBox (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
handles = updateImageListBox(handles);

% Update handles structure
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over imageStackListBox.
function imageStackListBox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to imageStackListBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4


% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on run and none of its controls.
function run_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over run.
function run_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
