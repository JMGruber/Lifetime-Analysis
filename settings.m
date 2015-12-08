function varargout = settings(varargin)
% SETTINGS MATLAB code for settings.fig
%      SETTINGS, by itself, creates a new SETTINGS or raises the existing
%      singleton*.
%
%      H = SETTINGS returns the handle to a new SETTINGS or the handle to
%      the existing singleton*.
%
%      SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SETTINGS.M with the given input arguments.
%
%      SETTINGS('Property','Value',...) creates a new SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help settings

% Last Modified by GUIDE v2.5 26-Jun-2014 15:05:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @settings_OpeningFcn, ...
                   'gui_OutputFcn',  @settings_OutputFcn, ...
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


% --- Executes just before settings is made visible.
function settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to settings (see VARARGIN)

handles.testQuit = 0;
if ~isempty(varargin)
	set(handles.txtCurDir,'String',varargin{1});
	set(handles.txtDirOut,'String',varargin{1});
else
	set(handles.txtCurDir,'String',pwd);
	set(handles.txtDirOut,'String',pwd);
	handles.testQuit = 1;
end

if exist([get(handles.txtCurDir,'String') '\settings.mat'],'file')
	load_settings(hObject,eventdata,handles);
	removeFileCheck(hObject, eventdata, handles);
end

% Choose default command line output for settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = [handles.output get(handles.txtCurDir,'String')];
varargout{1} = handles.output ;


% --- Executes on button press in btnNewDir.
function btnNewDir_Callback(hObject, eventdata, handles)
% hObject    handle to btnNewDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.txtCurDir,'String',uigetdir(get(handles.txtCurDir,'String'),'Choose a new directory'));

if exist([get(handles.txtCurDir,'String') '\settings.mat'],'file')
	load_settings(hObject,eventdata,handles);
end
removeFileCheck(hObject, eventdata, handles);



function txtFileName_Callback(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFileName as text
%        str2double(get(hObject,'String')) returns contents of txtFileName as a double


% --- Executes during object creation, after setting all properties.
function txtFileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFLast_Callback(hObject, eventdata, handles)
% hObject    handle to txtFLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFLast as text
%        str2double(get(hObject,'String')) returns contents of txtFLast as a double


% --- Executes during object creation, after setting all properties.
function txtFLast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFFirst_Callback(hObject, eventdata, handles)
% hObject    handle to txtFFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFFirst as text
%        str2double(get(hObject,'String')) returns contents of txtFFirst as a double


% --- Executes during object creation, after setting all properties.
function txtFFirst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtFSkip_Callback(hObject, eventdata, handles)
% hObject    handle to txtFSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtFSkip as text
%        str2double(get(hObject,'String')) returns contents of txtFSkip as a double


% --- Executes during object creation, after setting all properties.
function txtFSkip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtFSkip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnOutputFolder.
function btnOutputFolder_Callback(hObject, eventdata, handles)
% hObject    handle to btnOutputFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.txtDirOut,'String',uigetdir(get(handles.txtDirOut,'String'),'Choose a new output directory'));
removeFileCheck(hObject, eventdata, handles);

function txtChannel_Callback(hObject, eventdata, handles)
% hObject    handle to txtChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtChannel as text
%        str2double(get(hObject,'String')) returns contents of txtChannel as a double


% --- Executes during object creation, after setting all properties.
function txtChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chbHideFig.
function chbHideFig_Callback(hObject, eventdata, handles)
% hObject    handle to chbHideFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbHideFig


% --- Executes on button press in chbDrawLevels.
function chbDrawLevels_Callback(hObject, eventdata, handles)
% hObject    handle to chbDrawLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbDrawLevels


% --- Executes on button press in chbDensityMap.
function chbDensityMap_Callback(hObject, eventdata, handles)
% hObject    handle to chbDensityMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbDensityMap


% --- Executes on button press in chbParallel.
function chbParallel_Callback(hObject, eventdata, handles)
% hObject    handle to chbParallel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbParallel


% --- Executes on button press in chbForce.
function chbForce_Callback(hObject, eventdata, handles)
% hObject    handle to chbForce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbForce



function txtLiveThresh_Callback(hObject, eventdata, handles)
% hObject    handle to txtLiveThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLiveThresh as text
%        str2double(get(hObject,'String')) returns contents of txtLiveThresh as a double


% --- Executes during object creation, after setting all properties.
function txtLiveThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLiveThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtLiveTime_Callback(hObject, eventdata, handles)
% hObject    handle to txtLiveTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtLiveTime as text
%        str2double(get(hObject,'String')) returns contents of txtLiveTime as a double


% --- Executes during object creation, after setting all properties.
function txtLiveTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtLiveTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtThr2c_Callback(hObject, eventdata, handles)
% hObject    handle to txtThr2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtThr2c as text
%        str2double(get(hObject,'String')) returns contents of txtThr2c as a double


% --- Executes during object creation, after setting all properties.
function txtThr2c_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtThr2c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtIntBin_Callback(hObject, eventdata, handles)
% hObject    handle to txtIntBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIntBin as text
%        str2double(get(hObject,'String')) returns contents of txtIntBin as a double


% --- Executes during object creation, after setting all properties.
function txtIntBin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIntBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtTimeRes_Callback(hObject, eventdata, handles)
% hObject    handle to txtTimeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtTimeRes as text
%        str2double(get(hObject,'String')) returns contents of txtTimeRes as a double


% --- Executes during object creation, after setting all properties.
function txtTimeRes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtTimeRes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtGuessBg_Callback(hObject, eventdata, handles)
% hObject    handle to txtGuessBg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtGuessBg as text
%        str2double(get(hObject,'String')) returns contents of txtGuessBg as a double


% --- Executes during object creation, after setting all properties.
function txtGuessBg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtGuessBg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtIntFact_Callback(hObject, eventdata, handles)
% hObject    handle to txtIntFact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtIntFact as text
%        str2double(get(hObject,'String')) returns contents of txtIntFact as a double


% --- Executes during object creation, after setting all properties.
function txtIntFact_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtIntFact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtHistBin_Callback(hObject, eventdata, handles)
% hObject    handle to txtHistBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtHistBin as text
%        str2double(get(hObject,'String')) returns contents of txtHistBin as a double


% --- Executes during object creation, after setting all properties.
function txtHistBin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtHistBin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtQextent_Callback(hObject, eventdata, handles)
% hObject    handle to txtQextent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtQextent as text
%        str2double(get(hObject,'String')) returns contents of txtQextent as a double


% --- Executes during object creation, after setting all properties.
function txtQextent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtQextent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chbTestPhoton.
function chbTestPhoton_Callback(hObject, eventdata, handles)
% hObject    handle to chbTestPhoton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbTestPhoton



function txtNf_Callback(hObject, eventdata, handles)
% hObject    handle to txtNf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNf as text
%        str2double(get(hObject,'String')) returns contents of txtNf as a double


% --- Executes during object creation, after setting all properties.
function txtNf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtPhotonFactor_Callback(hObject, eventdata, handles)
% hObject    handle to txtPhotonFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtPhotonFactor as text
%        str2double(get(hObject,'String')) returns contents of txtPhotonFactor as a double


% --- Executes during object creation, after setting all properties.
function txtPhotonFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtPhotonFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chbFixBg.
function chbFixBg_Callback(hObject, eventdata, handles)
% hObject    handle to chbFixBg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbFixBg


% --- Executes on button press in chbTestSM.
function chbTestSM_Callback(hObject, eventdata, handles)
% hObject    handle to chbTestSM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbTestSM


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

temp = get(handles.txtFFirst,'String');
allfiles = str2double(get(handles.txtFFirst,'String')):str2double(get(handles.txtFLast,'String'));
temp = get(handles.txtFSkip,'String');
if ~isempty(temp)
spaces = strfind(temp,' ');
skipfiles = [];
if ~isempty(spaces)
	for i = length(spaces):-1:1
		skipfiles = [str2double(temp(spaces(i):end)) skipfiles];
		temp = temp(1:spaces(i)-1);
	end
end
skipfiles = [str2double(temp(1:end)) skipfiles];
skipfiles = sort(skipfiles);
else
	skipfiles = [];
end
pathname = get(handles.txtCurDir,'String');
readdir = pathname;
handles.output = pathname;
writedir = get(handles.txtDirOut,'String');
filename2 = get(handles.txtFileName,'String');
routerChannel = str2double(get(handles.txtChannel,'String')); % 1 for old APD, 2 for new MPD
drawlevels = get(handles.chbDrawLevels,'Value');  % draw intensity levels
hidelevel = get(handles.chbHideFig,'Value');	% Hides figure generation for quicker run
drawdensitymaps = get(handles.chbDensityMap,'Value'); % draw intensity maps
Parallel = get(handles.chbParallel,'Value');	% Then uses PARFOR
forceRun = get(handles.chbForce,'Value'); % To force code to try and skip all problematic traces and complete run. But runs twice.
Nf = str2double(get(handles.txtNf,'String')); % noise multiplication factor (1.5 - 2 for good S/N)
thrliveI = str2double(get(handles.txtLiveThresh,'String'));  % threshold for final intensity level (in cps, background included, intbin excluded)
thrlivet = str2double(get(handles.txtLiveTime,'String'));  % minimal survival time (in seconds)
thr2c = str2double(get(handles.txtThr2c,'String'));  % intensity threshold when definitely >1 complex (in cps, background included, intbin excluded)
intbin = str2double(get(handles.txtIntBin,'String'));   % number of intensity values to be binned together (input)---- for new version the definition is slightly different
timeres = str2double(get(handles.txtTimeRes,'String')); % time resolution (excl intbin)
fixbg = get(handles.chbFixBg,'Value'); % usually more accurate
guessbg = str2double(get(handles.txtGuessBg,'String')); % estimated background; used to test for SM (in cps). If fixbg=false, set rather too large than to small
testSM = get(handles.chbTestSM,'Value'); % test for single quantum unit by <=2-time step into quenched state
testphotonburst = get(handles.chbTestPhoton,'Value');
photonburstfactor = str2double(get(handles.txtPhotonFactor,'String'));   % factor by which dwelling time in Q state exceeds dwelling time in unQ states to define photonburst (typically 10)
intermedfactor = str2double(get(handles.txtIntFact,'String')); % parameter used to test if there are too many unnatural fluctuations (signifying unstable complex)
% decrease if more fluctuations should be incl. (e.g. 0.9)
histbin = str2double(get(handles.txtHistBin,'String'));  % number of intensity levels to be binned for plots (output)
Qextent = str2double(get(handles.txtQextent,'String'));  % extent of quenching (for panel 4)
highQ = get(handles.chbHQtraces,'Value'); % save high quality trace images
if length(get(handles.chbRemovePR,'Enable')) == 3
	removePR = 0;
else
	removePR = get(handles.chbRemovePR,'Value'); % remove previous results
end
save([readdir '\settings.mat'],'allfiles','skipfiles','pathname',...
	'readdir','writedir','filename2','routerChannel','drawlevels',...
	'hidelevel','drawdensitymaps','Parallel','forceRun','Nf','thrliveI',...
	'thrlivet','thr2c','intbin','timeres','fixbg','guessbg','testSM',...
	'testphotonburst','photonburstfactor','intermedfactor','histbin',...
	'Qextent','highQ','removePR');


% --- Executes on button press in btnContinue.
function btnContinue_Callback(hObject, eventdata, handles)
% hObject    handle to btnContinue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% uiresume(gcbf);
figure1_CloseRequestFcn(hObject, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(gcbf);

if handles.testQuit
	delete(gcbf);
end

function load_settings(hObject, eventdata, handles)
% My own function for loading directroy-saved settings.

load([get(handles.txtDirOut,'String') '\settings.mat']);

set(handles.txtFFirst,'String',int2str(allfiles(1)));
set(handles.txtFLast,'String',int2str(allfiles(end)));
temp = '';
skipfiles = sort(skipfiles);
for i = 1:length(skipfiles)
	if i == 1
		temp = int2str(skipfiles(i));
	else
		temp = [temp ' ' int2str(skipfiles(i))];
	end
end
set(handles.txtFSkip,'String',temp);
set(handles.txtCurDir,'String',pathname);
set(handles.txtDirOut,'String',writedir);
set(handles.txtFileName,'String', filename2);
set(handles.txtChannel,'String',int2str(routerChannel))
set(handles.chbDrawLevels,'Value',drawlevels);
set(handles.chbHideFig,'Value',hidelevel);
set(handles.chbDensityMap,'Value',drawdensitymaps);
set(handles.chbParallel,'Value',Parallel);
set(handles.chbForce,'Value',forceRun);
set(handles.txtNf,'String', num2str(Nf));
set(handles.txtLiveThresh,'String',num2str(thrliveI));
set(handles.txtLiveTime,'String',num2str(thrlivet));
set(handles.txtThr2c,'String',num2str(thr2c));
set(handles.txtIntBin,'String',num2str(intbin));
set(handles.txtTimeRes,'String',num2str(timeres));
set(handles.chbFixBg,'Value',fixbg);
set(handles.txtGuessBg,'String',num2str(guessbg));
set(handles.chbTestSM,'Value',testSM);
set(handles.chbTestPhoton,'Value',testphotonburst);
set(handles.txtPhotonFactor,'String',num2str(photonburstfactor));
set(handles.txtIntFact,'String',num2str(intermedfactor));
set(handles.txtHistBin,'String',num2str(histbin));
set(handles.txtQextent,'String',num2str(Qextent));
set(handles.chbHQtraces,'Value',highQ);
set(handles.chbRemovePR,'Value',removePR);


% --- Executes on button press in chbHQtraces.
function chbHQtraces_Callback(hObject, eventdata, handles)
% hObject    handle to chbHQtraces (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbHQtraces


% --- Executes on button press in chbRemovePR.
function chbRemovePR_Callback(hObject, eventdata, handles)
% hObject    handle to chbRemovePR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chbRemovePR

function removeFileCheck(hObject, eventdata, handles)
workingDir = get(handles.txtCurDir,'String');
resultDir = get(handles.txtDirOut,'String');
if length(workingDir) == length(resultDir)
	if workingDir == resultDir
		set(handles.chbRemovePR,'Enable','off');
	else
		set(handles.chbRemovePR,'Enable','on');
	end
else
	set(handles.chbRemovePR,'Enable','on');
end
