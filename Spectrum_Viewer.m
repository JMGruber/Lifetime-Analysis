function varargout = Spectrum_Viewer(varargin)
% SPECTRUM_VIEWER MATLAB code for Spectrum_Viewer.fig
%      SPECTRUM_VIEWER, by itself, creates a new SPECTRUM_VIEWER or raises the existing
%      singleton*.
%
%      H = SPECTRUM_VIEWER returns the handle to a new SPECTRUM_VIEWER or the handle to
%      the existing singleton*.
%
%      SPECTRUM_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTRUM_VIEWER.M with the given input arguments.
%
%      SPECTRUM_VIEWER('Property','Value',...) creates a new SPECTRUM_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spectrum_Viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spectrum_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spectrum_Viewer

% Last Modified by GUIDE v2.5 22-Sep-2014 10:37:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Spectrum_Viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @Spectrum_Viewer_OutputFcn, ...
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




% --- Executes just before Spectrum_Viewer is made visible.
function Spectrum_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spectrum_Viewer (see VARARGIN)

% Choose default command line output for Spectrum_Viewer
handles.output = hObject;
handles.plot=hObject;
handles.folder_name=hObject;
handles.Meanval=hObject;
handles.Meanval=[];
handles.YLim=hObject;


% Update handles structure
set(handles.text2, 'String', 'Directory');
set(handles.edit2, 'String', '200');
set(handles.edit5, 'String', '8000');
set(handles.slider2,'Min',0);
set(handles.slider2,'Max',200);
set(handles.slider2,'Value',1);
set(handles.slider2,'SliderStep',[1/(200-1) 1/(200-1)]);
set(handles.checkbox1,'Value',1);
set(handles.slider4,'Value',1);
set(handles.slider4,'Min',1);
set(handles.slider4,'Max',60);
set(handles.checkbox2,'Value',0);
handles.folder_name='O:\Michael\2014';

guidata(hObject, handles);






% UIWAIT makes Spectrum_Viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Spectrum_Viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.text2;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
j=get(handles.slider2,'Value');
handles.folder_name=uigetdir(strcat(handles.folder_name), 'Pick a Directory');
guidata(hObject, handles);
index_selected = get(handles.listbox3,'Value');
list = get(handles.listbox3,'String');
item_selected = list{index_selected};
set(handles.text2, 'String', strcat(handles.folder_name,'\', item_selected));

errorcheck(j,handles.text2);
guidata(hObject, handles);
slider2_Callback(hObject, eventdata, handles); 

%guidata(hObject, handles);



function errorcheck(nummer,dir)
if ~exist(strcat(get(dir,'String'),int2str(nummer)));    
    cla;
    display('No such Spectrum')
    error('No such Spectrum')
    
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


cla;
if get(handles.checkbox5,'Value')==0
subplot(1,1,1);
set(gca, 'Position',[0.03 0.15  0.85 0.65])

end

j=get(handles.slider2,'Value');
k=round(get(handles.slider4,'Value'));
set(handles.text3,'String',strcat('File Nr. ',int2str(j)));

    if strcmp(get(handles.text2,'String'),'Directory')||strcmp(get(handles.text2,'String'),'0');
        pushbutton1_Callback(hObject, eventdata, handles)
    end
    
%     if exist(strcat(get(handles.text2,'String'),int2str(j)),'file')==0
%         pushbutton1_Callback(hObject, eventdata, handles)
%     end
    
    
errorcheck(j,handles.text2);
mat=dlmread(strcat(get(handles.text2,'String'),int2str(j)));

 file_list = get(handles.listbox3,'String');
 index_selected = get(handles.listbox3,'Value');
 test = file_list{index_selected};

if strcmp(test,'Decaytrace'),
else
if get(handles.checkbox3,'Value')==1&&~isempty(handles.Meanval);
 for i=2:size(mat,2)
     mat(:,i)=mat(:,i)-handles.Meanval;
 end
end
end




    if get(handles.checkbox6,'Value')==1        
    subplot(2,2,1)
    set(gca, 'Position',[0.05 0.45  0.40 0.3])
    end 
    
    if get(handles.checkbox5,'Value')==1
        
    subplot(2,2,1)
    set(gca, 'Position',[0.05 0.45  0.40 0.3])


    end    

    if get(handles.checkbox1,'Value')==1
        handles.plot=plot(mat(:,1),mat(:,2:length(mat(1,:))));
            if get(handles.checkbox2,'Value')==1
        
            set(gca, 'YLim',[handles.YLim]);
            

            end
        handles.YLim=get(gca,'YLim');
        

    else
        handles.plot=plot(mat(:,1),mat(:,1+k)); 
            if get(handles.checkbox2,'Value')==1
        
            set(gca, 'YLim',[handles.YLim]);
            end
        handles.YLim=get(gca,'YLim');


    end


if strcmp(test,'Decaytrace'),
    set(gca,'YScale','log');
end
    
   if get(handles.checkbox5,'Value')==1 


subplot('Position',[0.50 0.1  0.40 0.3]);
mat2=dlmread(strcat(handles.folder_name,'\trace',int2str(j)));
plot(mat2(:,1),mat2(:,2));

ylabel('Fluorescence [cps]') ;
xlabel('Time [s]');
b=str2double(get(handles.edit5,'String'));
ylim([0 b]);
subplot('Position',[0.05 0.1  0.40 0.3]);
mat2=dlmread(strcat(handles.folder_name,'\decaytrace',int2str(j)));
plot(mat2(:,1),mat2(:,2));
ylabel('Fluorescence [cps]');
xlabel('Time [s]');
set(gca,'YScale','log'); %test
   end
   
      if get(handles.checkbox6,'Value')==1 
        subplot('Position',[0.50 0.45  0.40 0.3]);
        image(mat(:,1),[1:(length(mat(1,:))-1)],mat(:,2:length(mat(1,:)))','CDataMapping','Scaled')%,'CDataMapping','scaled'
        set(gcf,'name',strcat('Spectrum Nr. ',int2str(j)));
        set(gca,'YDir','normal');
        caxis([-20 150]);
        ylabel('Time [s]') ;
        xlabel('Wavelength [nm]');
      end
      if get(handles.checkbox7,'Value')==1 
        subplot('Position',[0.05 0.1  0.40 0.3]);
            for i=2:length(mat(1,:))
                d(i)=sum(mat(1:length(mat(:,1)),i));
            end
        plot(d(2:length(mat(1,:))));
        set(gcf,'name',strcat('Spectrum Nr. ',int2str(j)));
        set(gca,'YDir','normal');
        ylabel('Integrated Intensity [cps]') ;
        xlabel('Time [s]');
      end
  


    
hold off
guidata(hObject, handles);
if get(handles.checkbox4,'Value')==1
   pushbutton7_Callback(hObject, eventdata, handles);
    
end





% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
a=str2double(get(hObject,'String'));
set(handles.slider2,'Value',1);
set(handles.slider2,'Max',a);
set(handles.slider2,'SliderStep',[1/(a-1) 1/(a-1)]);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if ~strcmp(get(handles.text2,'String'),'Directory')||strcmp(get(handles.text2,'String'),'0');
    if get(hObject,'Value')==1
        slider2_Callback(hObject, eventdata, handles);
    else
        slider4_Callback(hObject, eventdata, handles)
    end
end



% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if strcmp(get(handles.text2,'String'),'Directory')||strcmp(get(handles.text2,'String'),'0');
    pushbutton1_Callback(hObject, eventdata, handles)
end

set(handles.checkbox1,'Value',0);
j=get(handles.slider2,'Value');
errorcheck(j,handles.text2);
mat=dlmread(strcat(get(handles.text2,'String'),int2str(j)));
maximum=size(mat,2)-1;
set(handles.slider4,'Min',1);
    if maximum < get(handles.slider4,'Value')
       set(handles.slider4,'Value',maximum);
    end
set(handles.slider4,'Max',maximum);
set(handles.slider4,'SliderStep',[1/(maximum-1) 1/(maximum-1)]);
i=round(get(handles.slider4,'Value'));

if get(handles.checkbox3,'Value')==1&&~isempty(handles.Meanval);
for k=2:size(mat,2)
    mat(:,k)=mat(:,k)-handles.Meanval;
end
end

if get(handles.checkbox2,'Value')==1
    cla;
   hold on
    axis manual;
end
    
handles.plot=plot(mat(:,1),mat(:,i+1));
hold off
set(handles.text5,'String',strcat('Spectrum Nr. ',int2str(i)));
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2

slider2_Callback(hObject, eventdata, handles);







% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
j=get(handles.slider2,'Value');
dir=get(handles.text2,'String');


 if get(handles.checkbox4,'Value')==1
     figure(1);

 else
     figure();
 end
mat=dlmread(strcat(dir,int2str(j)));
    if get(handles.checkbox3,'Value')==1

    
for k=2:size(mat,2)
    mat(:,k)=mat(:,k)-handles.Meanval;
end
    end

       if get(handles.checkbox6,'Value')==1 

        image(mat(:,1),[1:(length(mat(1,:))-1)],mat(:,2:length(mat(1,:)))','CDataMapping','Scaled')%,'CDataMapping','scaled'
        set(gcf,'name',strcat('Spectrum Nr. ',int2str(j)));
        set(gca,'YDir','normal');
        %caxis([-20 200]);
        ylabel('Time [s]') ;
        xlabel('Wavelength [nm]');
      end
      if get(handles.checkbox7,'Value')==1 

            for i=2:length(mat(1,:))
                d(i)=sum(mat(1:length(mat(:,1)),i));
            end
        plot(d(2:length(mat(1,:))));
        set(gcf,'name',strcat('Spectrum Nr. ',int2str(j)));
        set(gca,'YDir','normal');
        ylabel('Integrated Intensity [cps]') ;
        xlabel('Time [s]');
      end   
    
    
    
    
    
    
    
% image(mat(:,1),[1:(length(mat(1,:))-1)],mat(:,2:length(mat(1,:)))','CDataMapping','scaled')
% set(gcf,'name',strcat('Spectrum Nr. ',int2str(j)));
% set(gca,'YDir','normal');
% ylabel('Time [s]') ;
% xlabel('Wavelength [nm]');






% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(3);
hold all

Xaxis=get(handles.plot,'Xdata');
Yaxis=get(handles.plot,'ydata');

if iscell(Xaxis)
    Xmat=[];
    Ymat=[];
    for i=1:length(Xaxis)
        Xmat(:,i)=Xaxis{i,1}';
        Ymat(:,i)=Yaxis{i,1}';
    end
    plot(Xmat(:,1:length(Xaxis)),Ymat(:,1:length(Yaxis)));
    ylabel('Fluorescence [cps]') ;
xlabel('Wavelength [nm]');
    
else
    
     Ymax=max(Yaxis);
     %Yaxis=(Yaxis)./Ymax+0.01;
    plot(Xaxis,Yaxis);
    ylabel('Fluorescence [cps]') ;
xlabel('Wavelength [nm]');
end
hold off

 file_list = get(handles.listbox3,'String');
 index_selected = get(handles.listbox3,'Value');
 test = file_list{index_selected};
if strcmp(test,'Decaytrace'),
    set(gca,'YScale','log');
    ylabel('Counts') ;
    xlabel('Delay time[ns]');
    axis([0 12 1 1000])
end

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
j=get(handles.slider2,'Value');
mat=dlmread(strcat(get(handles.text2,'String'),int2str(j)));
%Meanval=dlmread('D:\Sheela\Mutant 3 5C POC flushed e-5\bg.txt');
Meanval=mean(mat(:,2:size(mat,2)),2);
    if get(handles.checkbox2,'Value')==1
        cla;
       hold on
        axis manual;
    end
handles.Meanval=Meanval;
handles.plot=plot(mat(:,1),Meanval);
hold off
guidata(hObject, handles);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
slider2_Callback(hObject, eventdata, handles);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, ~, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure();

j=get(handles.slider2,'Value');
mat=dlmread(strcat(handles.folder_name,'\trace',int2str(j)));
plot(mat(:,1),mat(:,2));
ylabel('Fluorescence [cps]') ;
xlabel('Time [s]');




% --- Executes on selection change in listbox3.
function listbox3_Callback(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox3
index_selected = get(handles.listbox3,'Value');
list = get(handles.listbox3,'String');
item_selected = list{index_selected};

set(handles.text2, 'String', strcat(handles.folder_name,'\', item_selected));
slider2_Callback(hObject, eventdata, handles); 




% --- Executes during object creation, after setting all properties.
function listbox3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
