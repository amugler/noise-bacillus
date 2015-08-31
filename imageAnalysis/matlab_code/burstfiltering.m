function varargout = burstfiltering(varargin)
% BURSTFILTERING M-file for burstfiltering.fig
%      BURSTFILTERING, by itself, creates a new BURSTFILTERING or raises the existing
%      singleton*.
%
%      BURSTFILTERING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BURSTFILTERING.M with the given input arguments.
%
%      BURSTFILTERING('Property','Value',...) creates a new BURSTFILTERING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before burstfiltering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to burstfiltering_OpeningFcn via varargin.
%
%      INPUT
%      image : Image data 2D matrix, double
%
%      Parameters
%      cc     : connected components structure (see bwconncomp)
%      rp     : region properties (see regionprops)
%      ranges : structure specifying ranges (1st output of this function)
%      mask   : binary mask to analyze regions marked as true or 1
%
%      OUTPUT
%      ranges : structure specifying the property ranges selected
%      cc     : connected components structure
%      rp     : region properties structure
%      output : handle to the figure
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help burstfiltering

% Last Modified by GUIDE v2.5 20-Oct-2010 14:30:22

% Mark Kittisopikul, 2010
% UT Southwestern

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @burstfiltering_OpeningFcn, ...
                   'gui_OutputFcn',  @burstfiltering_OutputFcn, ...
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

% --- Outputs from this function are returned to the command line.
function varargout = burstfiltering_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.ranges;
varargout{2} = handles.burstCC;
varargout{3} = handles.burstRP;
varargout{4} = handles.output;
close(handles.figure1);

% --- Executes just before burstfiltering is made visible.
function burstfiltering_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to burstfiltering (see VARARGIN)

% Choose default command line output for burstfiltering
handles.output = hObject;

parser = inputParser;
parser.addOptional('image',[]);
% parser.addParamValue('filter',[]);
parser.addParamValue('cc',[]);
parser.addParamValue('rp',[]);
parser.addParamValue('ranges',[]);
parser.addParamValue('mask',[]);
parser.parse(varargin{:});
I = parser.Results.image;
if(isempty(parser.Results.cc))
    CC = edgeDetect(I);
else
    CC = parser.Results.cc;
end
orgLabels = labelmatrix(CC);
if(isempty(parser.Results.rp))
    RP = getRegionProps(I,CC);
else
    RP = parser.Results.rp;
end

labelColorMap = jet(CC.NumObjects);
[s, idx] = sort(rand(1,CC.NumObjects));
labelColorMap = labelColorMap(idx,:);
labelColorMap = [ [0 0 0] ; labelColorMap ];

if(isempty(parser.Results.ranges))
    handles.ranges = getRanges(RP);
else
    handles.ranges = parser.Results.ranges;
end
handles.sorted = getSorted(RP);

handles.maskFilter = ones(1,CC.NumObjects);
if(isempty(parser.Results.mask))
else
    outOfMask = unique(orgLabels(~parser.Results.mask));
    handles.maskFilter(outOfMask(outOfMask ~= 0)) = 0;
end

propertyName = get(get(handles.uipanel1,'SelectedObject'),'UserData');
maxValue = length(handles.sorted.(propertyName));

set(handles.minSlider,'Max',maxValue-1);
set(handles.minSlider,'Min',1);
set(handles.minSlider,'Value',1);
set(handles.minSlider,'SliderStep',[1/(maxValue-1) 0.1]);
set(handles.maxSlider,'Max',maxValue);
set(handles.maxSlider,'Min',2);
set(handles.maxSlider,'SliderStep',[1/(maxValue-1) 0.1]);
set(handles.maxSlider,'Value',maxValue);

imshow(orgLabels,[],'colormap',labelColorMap);


updateUIRanges(handles,handles.ranges);

propNameMap =struct();
propNameMap.area = 'Area';
propNameMap.min = 'MinIntensity';
propNameMap.mean = 'MeanIntensity';
propNameMap.max = 'MaxIntensity';
propNameMap.solidity = 'Solidity';
propNameMap.euler = 'EulerNumber';
propNameMap.extent = 'Extent';
propNameMap.perimeter = 'Perimeter';
propNameMap.eccentricity = 'Eccentricity';

handles.propNameMap = propNameMap;
handles.CC = CC;
handles.RP = RP;
handles.burstCC = CC;
handles.burstRP = RP;
handles.labelColorMap = labelColorMap;
handles.orgLabels = orgLabels;
handles.burstLabels = orgLabels;
handles.propertyName = propertyName;
handles.I = I;

% Update handles structure
guidata(hObject, handles);

if(~isempty(parser.Results.ranges))
    updateFilter(hObject,handles);
    handles = guidata(hObject);
%     viewBursts(handles);
    updateSliders(handles);
end

% UIWAIT makes burstfiltering wait for user response (see UIRESUME)
uiwait(handles.figure1);

 

    
function [RP] = getRegionProps(I,CC) 
% regionproperties
RP = regionprops(CC,I,'MeanIntensity','MinIntensity','MaxIntensity','Area','Solidity','Extent','EulerNumber','Perimeter','Eccentricity');

function updateUIRanges(handles,ranges)
set(handles.areaMin,'String',ranges.area(1));
set(handles.areaMax,'String',ranges.area(2));
set(handles.minMin,'String',ranges.min(1));
set(handles.minMax,'String',ranges.min(2))
set(handles.meanMin,'String',ranges.mean(1));
set(handles.meanMax,'String',ranges.mean(2));
set(handles.maxMin,'String',ranges.max(1));
set(handles.maxMax,'String',ranges.max(2));
set(handles.solidityMin,'String',ranges.solidity(1));
set(handles.solidityMax,'String',ranges.solidity(2));
set(handles.extentMin,'String',ranges.extent(1));
set(handles.extentMax,'String',ranges.extent(2));
set(handles.eulerMin,'String',ranges.euler(1));
set(handles.eulerMax,'String',ranges.euler(2));
set(handles.perimeterMin,'String',ranges.perimeter(1));
set(handles.perimeterMax,'String',ranges.perimeter(2));
set(handles.eccentricityMin,'String',ranges.eccentricity(1));
set(handles.eccentricityMax,'String',ranges.eccentricity(2));

function ranges = getUIRanges(handles)
ranges = struct();
ranges.area(1) = str2double(get(handles.areaMin,'String'));
ranges.area(2) = str2double(get(handles.areaMax,'String'));
ranges.min(1) = str2double(get(handles.minMin,'String'));
ranges.min(2) = str2double(get(handles.minMax,'String'));
ranges.mean(1) = str2double(get(handles.meanMin,'String'));
ranges.mean(2) = str2double(get(handles.meanMax,'String'));
ranges.max(1) = str2double(get(handles.maxMin,'String'));
ranges.max(2) = str2double(get(handles.maxMax,'String'));
ranges.solidity(1) = str2double(get(handles.solidityMin,'String'));
ranges.solidity(2) = str2double(get(handles.solidityMax,'String'));
ranges.extent(1) = str2double(get(handles.extentMin,'String'));
ranges.extent(2) = str2double(get(handles.extentMax,'String'));
ranges.euler(1) = str2double(get(handles.eulerMin,'String'));
ranges.euler(2) = str2double(get(handles.eulerMax,'String'));
ranges.perimeter(1) = str2double(get(handles.perimeterMin,'String'));
ranges.perimeter(2) = str2double(get(handles.perimeterMax,'String'));
ranges.eccentricity(1) = str2double(get(handles.eccentricityMin,'String'));
ranges.eccentricity(2) = str2double(get(handles.eccentricityMax,'String'));


function [range] =getRange(v)
range = [floor(min(v)) ceil(max(v))];


function ranges = getRanges(RP)
ranges.area = getRange([RP.Area]);
ranges.min = getRange([RP.MinIntensity]);
ranges.mean = getRange([RP.MeanIntensity]);
ranges.max = getRange([RP.MaxIntensity]);
ranges.solidity = getRange([RP.Solidity]);
ranges.extent = getRange([RP.Extent]);
ranges.euler = getRange([RP.EulerNumber]);
ranges.perimeter = getRange([RP.Perimeter]);
ranges.eccentricity = getRange([RP.Eccentricity]);

function values = getSorted(RP)
values.area = [-Inf unique(sort([RP.Area])) Inf];
values.min = [-Inf unique(sort([RP.MinIntensity])) Inf];
values.mean = [-Inf unique(sort([RP.MeanIntensity])) Inf];
values.max = [-Inf unique(sort([RP.MaxIntensity])) Inf];
values.solidity = [-Inf unique(sort([RP.Solidity])) Inf];
values.extent = [-Inf unique(sort([RP.Extent])) Inf];
values.euler = [-Inf unique(sort([RP.EulerNumber])) Inf];
values.perimeter = [-Inf unique(sort([RP.Perimeter])) Inf];
values.eccentricity = [-Inf unique(sort([RP.Eccentricity])) Inf];

function updateFilter(hObject,handles)
%     disp('updating filter');
ranges = handles.ranges;
RP = handles.RP;
CC = handles.CC;
% burstRP = handles.burstRP;
% burstCC = handles.burstCC;
F = zeros(10,CC.NumObjects);
F(1,:) = rangeFilter([RP.Area],ranges.area);
F(2,:) = rangeFilter([RP.MinIntensity],ranges.min);
F(3,:) = rangeFilter([RP.MeanIntensity],ranges.mean);
F(4,:) = rangeFilter([RP.MaxIntensity],ranges.max);
F(5,:) = rangeFilter([RP.Solidity],ranges.solidity);
F(6,:) = rangeFilter([RP.Extent],ranges.extent);
F(7,:) = rangeFilter([RP.EulerNumber],ranges.euler);
F(8,:) = rangeFilter([RP.Perimeter],ranges.perimeter);
F(9,:) = rangeFilter([RP.Eccentricity],ranges.eccentricity);
F(10,:) = handles.maskFilter;

filter = all(F,1);
burstCC = filtercc(CC,filter);
burstRP = RP(filter);
labelMatrix = labelmatrix(burstCC);
% imshow(labelMatrix,[],'colormap',handles.labelColorMap);


handles.burstCC = burstCC;
handles.burstRP = burstRP;
handles.burstLabels = labelMatrix;
viewBursts(handles);
guidata(hObject,handles);





% --- Executes on slider movement.
function minSlider_Callback(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minValue=round(get(hObject,'Value'));
maxValue=round(get(handles.maxSlider,'Value'));
if(minValue > maxValue)
    maxValue = minValue +1;
    set(handles.maxSlider,'Value',maxValue);
end
pN = handles.propertyName;
sorted = handles.sorted.(pN);
set(handles.([pN 'Min']),'String',sorted(minValue));
set(handles.([pN 'Max']),'String',sorted(maxValue));
handles.ranges = getUIRanges(handles);
guidata(hObject,handles);
updateFilter(hObject,handles);



% --- Executes during object creation, after setting all properties.
function minSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function maxSlider_Callback(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
maxValue=round(get(hObject,'Value'));
minValue=round(get(handles.minSlider,'Value'));
if(minValue > maxValue)
    minValue = maxValue -1;
    set(handles.minSlider,'Value',minValue);
end
pN = handles.propertyName;
sorted = handles.sorted.(pN);
set(handles.([pN 'Min']),'String',sorted(minValue));
set(handles.([pN 'Max']),'String',sorted(maxValue));
handles.ranges = getUIRanges(handles);
guidata(hObject,handles);
updateFilter(hObject,handles);

% --- Executes during object creation, after setting all properties.
function maxSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in viewMenu.
function viewMenu_Callback(hObject, eventdata, handles)
% hObject    handle to viewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns viewMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from viewMenu
viewBursts(handles);

function viewBursts(handles)
% disp('viewBursts');
contents = get(handles.viewMenu,'String');
selected =  contents{get(handles.viewMenu,'Value')};
burstrp = handles.burstRP;
name = handles.propNameMap.(handles.propertyName);
switch selected
    case 'Selected Bursts'
        imshow(handles.burstLabels,[],'colormap',handles.labelColorMap);
    case 'Outlines'
        imshow(makergb(edge(handles.burstLabels,'log',0),handles.I))
    case 'Property'
        burstrp = handles.burstRP;
        name = handles.propNameMap.(handles.propertyName);
        propmap = [burstrp.(name)];
        propmap = round((propmap - min(propmap))/(max(propmap) - min(propmap)) * 254) +2;
        propmap = [1 propmap];
        propimg = propmap(handles.burstLabels+1);
        imshow(propimg,[],'colormap',[[0 0 0];jet(255)]);
    case 'Histogram'
        hist(double([burstrp.(name)]));
end

% --- Executes during object creation, after setting all properties.
function viewMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to viewMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in identify.
function identify_Callback(hObject, eventdata, handles)
% hObject    handle to identify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [x,y] = ginput(1)
    x = round(x);
    y = round(y);
%     h = impoint;
    RP = handles.RP;
    orgLabels = handles.orgLabels;
    % Update position in title using newPositionCallback
%     addNewPositionCallback(h,@(h) set(gcf,'Name',sprintf('(%1.0f,%1.0f,%1.0d)',h(1),h(2),orgLabels(h(2),h(1)))));
%     addNewPositionCallback(h,@(h) disp(RP(orgLabels(h(2),h(1)))));
%     % Construct boundary constraint function
%     fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
%     % Enforce boundary constraint function using setPositionConstraintFcn
%     setPositionConstraintFcn(h,fcn);
%     setColor(h,'r');
%     p = round(wait(h));
%     orgIdx = orgLabels(p(2),p(1));
    orgIdx = orgLabels(y,x);
    if(orgIdx ~= 0)
        disp(RP(orgIdx));
%         keyboard
    end
%     delete(h);

function textfieldUpdate(hObject, eventdata, handles)
%asdf
handles.ranges = getUIRanges(handles);
guidata(hObject,handles);
updateFilter(hObject,handles);
updateSliders(handles);

function areaMin_Callback(hObject, eventdata, handles)
% hObject    handle to areaMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of areaMin as text
%        str2double(get(hObject,'String')) returns contents of areaMin as a double
textfieldUpdate(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function areaMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to areaMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function areaMax_Callback(hObject, eventdata, handles)
% hObject    handle to areaMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of areaMax as text
%        str2double(get(hObject,'String')) returns contents of areaMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function areaMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to areaMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minMin_Callback(hObject, eventdata, handles)
% hObject    handle to minMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minMin as text
%        str2double(get(hObject,'String')) returns contents of minMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function minMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minMax_Callback(hObject, eventdata, handles)
% hObject    handle to minMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minMax as text
%        str2double(get(hObject,'String')) returns contents of minMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function minMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function meanMin_Callback(hObject, eventdata, handles)
% hObject    handle to meanMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meanMin as text
%        str2double(get(hObject,'String')) returns contents of meanMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function meanMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function meanMax_Callback(hObject, eventdata, handles)
% hObject    handle to meanMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of meanMax as text
%        str2double(get(hObject,'String')) returns contents of meanMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function meanMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to meanMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxMin_Callback(hObject, eventdata, handles)
% hObject    handle to maxMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxMin as text
%        str2double(get(hObject,'String')) returns contents of maxMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function maxMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxMax_Callback(hObject, eventdata, handles)
% hObject    handle to maxMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxMax as text
%        str2double(get(hObject,'String')) returns contents of maxMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function maxMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function solidityMin_Callback(hObject, eventdata, handles)
% hObject    handle to solidityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of solidityMin as text
%        str2double(get(hObject,'String')) returns contents of solidityMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function solidityMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solidityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function solidityMax_Callback(hObject, eventdata, handles)
% hObject    handle to solidityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of solidityMax as text
%        str2double(get(hObject,'String')) returns contents of solidityMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function solidityMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solidityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eulerMin_Callback(hObject, eventdata, handles)
% hObject    handle to eulerMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eulerMin as text
%        str2double(get(hObject,'String')) returns contents of eulerMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function eulerMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eulerMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eulerMax_Callback(hObject, eventdata, handles)
% hObject    handle to eulerMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eulerMax as text
%        str2double(get(hObject,'String')) returns contents of eulerMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function eulerMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eulerMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function extentMin_Callback(hObject, eventdata, handles)
% hObject    handle to extentMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extentMin as text
%        str2double(get(hObject,'String')) returns contents of extentMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function extentMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extentMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function extentMax_Callback(hObject, eventdata, handles)
% hObject    handle to extentMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of extentMax as text
%        str2double(get(hObject,'String')) returns contents of extentMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function extentMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to extentMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function perimeterMin_Callback(hObject, eventdata, handles)
% hObject    handle to perimeterMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of perimeterMin as text
%        str2double(get(hObject,'String')) returns contents of perimeterMin as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function perimeterMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perimeterMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function perimeterMax_Callback(hObject, eventdata, handles)
% hObject    handle to perimeterMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of perimeterMax as text
%        str2double(get(hObject,'String')) returns contents of perimeterMax as a double
textfieldUpdate(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function perimeterMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to perimeterMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pN = handles.propNameMap.(handles.propertyName);
hist([handles.burstRP.(pN)]);

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
propertyName = get(eventdata.NewValue,'UserData');
handles.propertyName = propertyName;
guidata(hObject,handles);
updateSliders(handles);
viewBursts(handles);

function updateSliders(handles)
propertyName = handles.propertyName;
range = handles.ranges.(propertyName);
sorted = handles.sorted.(propertyName);
[propMin,propMax] = getSliderValues(sorted,range);
maxValue = length(sorted);
set(handles.minSlider,'Max',maxValue-1,'Value',propMin,'SliderStep',[1/(maxValue-1) 0.1]);
% set(handles.minSlider,'Value',range(1));
set(handles.maxSlider,'Max',maxValue,'Value',propMax,'SliderStep',[1/(maxValue-1) 0.1]);
% propMin
% propMax
% maxValue
% set(handles.maxSlider,'Value',range(2));
% handles.propertyName = propertyName;
% guidata(hObject,handles);
% disp('update sliders');
% viewBursts(handles);

function [propMin,propMax] = getSliderValues(sorted,range)
propMin = find(sorted <= range(1),1,'last');
propMax = find(sorted >= range(2),1,'first');
if(isempty(propMin))
    propMin = 1;
end
if(isempty(propMax))
    propMax = length(sorted);
end



% --- Executes on button press in areaRadio.
function areaRadio_Callback(hObject, eventdata, handles)
% hObject    handle to areaRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of areaRadio


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imshow(handles.I,[]);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imshow(handles.orgLabels,[],'colormap',handles.labelColorMap);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imshow(makergb(edge(handles.orgLabels,'log',0),handles.I))



function eccentricityMin_Callback(hObject, eventdata, handles)
% hObject    handle to eccentricityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eccentricityMin as text
%        str2double(get(hObject,'String')) returns contents of eccentricityMin as a double


% --- Executes during object creation, after setting all properties.
function eccentricityMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eccentricityMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function eccentricityMax_Callback(hObject, eventdata, handles)
% hObject    handle to eccentricityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eccentricityMax as text
%        str2double(get(hObject,'String')) returns contents of eccentricityMax as a double


% --- Executes during object creation, after setting all properties.
function eccentricityMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eccentricityMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
% close(handles.figure1);
