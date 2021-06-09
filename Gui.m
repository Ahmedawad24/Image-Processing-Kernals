function varargout = Gui(varargin)
% GUI MATLAB code for Gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gui

% Last Modified by GUIDE v2.5 08-Jan-2021 19:30:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_OutputFcn, ...
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


% --- Executes just before Gui is made visible.
function Gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gui (see VARARGIN)

% Choose default command line output for Gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global a;
global a2;
global a3;

a=imread('Input.jpg');
a3=a;
a = rgb2gray(a);
a2=a;
axes(handles.axes1);
imshow(a);
axes(handles.axes2);
imshow(a2);
axes(handles.axes3);
imshow(a3);



% UIWAIT makes Gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in Gaborpushbutton.
function Gaborpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Gaborpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
in=a2;
image_resize=imresize(in, [160 160]);
image_resize=im2double(image_resize);

gamma=0.3; %aspect ratio
          psi=0; %phase
          i=2;
          s=0;
          index=-1;
          
           %orientation
           
          bw=2.8; %bandwidth or effective width
          lambda=3.5; % wavelength
          pi=180;
          Max = zeros(160,160);
          Thetas = zeros(160,160);
          ConF = zeros(160,160);
       for x=1:160
              for y=1:160
                  max=-999999;
                  s=0;
                  for theta=360:-45:0 
                      x_theta=image_resize(x,y)*cos(theta)+image_resize(x,y)*sin(theta);
                      y_theta=-image_resize(x,y)*sin(theta)+image_resize(x,y)*cos(theta);
                      gb(x,y)= exp(-(x_theta.^2/2*bw^2+ gamma^2*y_theta.^2/2*bw^2))*cos(2*pi/lambda*x_theta+psi);
                        if gb(x,y) > max
                            max = gb(x,y);
                            index=s;
                        end
                        s=s+1;
                  end
                  Max(x,y) = gb(x,y);
                  Thetas(x,y)=index*45;
              end
       end
          axes(handles.axes2);
          imshow(Max);
          a2= Max;
setappdata(0,'Max',Max);
% --- Executes on button press in ConFpushbutton.
function ConFpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ConFpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
in=a2;
image_resize=imresize(in, [160 160]);
image_resize=im2double(image_resize);

          gamma=0.3; %aspect ratio
          psi=0; %phase
          i=2;
          s=0;
          index=-1;
          
           %orientation
           
          bw=2.8; %bandwidth or effective width
          lambda=3.5; % wavelength
          pi=180;
          Max = zeros(160,160);
          Thetas = zeros(160,160);
          ConF = zeros(160,160);
          for x=1:160
              for y=1:160
                  max=-999999;
                  s=0;
                  for theta=360:-45:0 
                      x_theta=image_resize(x,y)*cos(theta)+image_resize(x,y)*sin(theta);
                      y_theta=-image_resize(x,y)*sin(theta)+image_resize(x,y)*cos(theta);
                      gb(x,y)= exp(-(x_theta.^2/2*bw^2+ gamma^2*y_theta.^2/2*bw^2))*cos(2*pi/lambda*x_theta+psi);
                        if gb(x,y) > max
                            max = gb(x,y);
                            index=s;
                        end
                        s=s+1;
                  end
                  Max(x,y) = gb(x,y);
                  Thetas(x,y)=index*45;
              end
          end
       
          Distance=0;
          for x=1:160
              for y=1:160
                  tot=0;
                  for theta=360:-45:0
                      x_theta=image_resize(x,y)*cos(theta)+image_resize(x,y)*sin(theta);
                      y_theta=-image_resize(x,y)*sin(theta)+image_resize(x,y)*cos(theta);
                      gb(x,y)= exp(-(x_theta.^2/2*bw^2+ gamma^2*y_theta.^2/2*bw^2))*cos(2*pi/lambda*x_theta+psi);
                      h = theta-Thetas(x,y);
                      if h < 0 
                          h = h*-1;
                      end
                      ConF(x,y) = ConF(x,y) + ( (theta-Thetas(x,y)) *(gb(x,y)-Max(x,y)).^2).^0.5;
                      
                  end
              end
          end
          axes(handles.axes2);
          imshow(ConF);
          a2=ConF;
          setappdata(0,'ConF',ConF);       
% --- Executes on button press in Resetpushbutton.
function Resetpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Resetpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
global a2;
global a3;
a=imread('Input.jpg');
a3=a;
a = rgb2gray(a);
a2=a;
axes(handles.axes1);
imshow(a);
axes(handles.axes2);
imshow(a2);
axes(handles.axes3);
imshow(a3);
% --- Executes on button press in Guassianpushbutton.
function Guassianpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Guassianpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
I = a2;
image_resize=imresize(I, [160 160]);
image_resize=im2double(image_resize);
N=5; %must be odd
sigma=1;
x=1:N;
X=exp(-(x-((N+1)/2)).^2/(2*sigma^2));
h=X'*X;
h=h./sum(h(:));
%I=filter2(h,I); %this is faster
[is,js]=size(image_resize);
Ib = NaN(is+N-1,js+N-1); %add borders
b=(N-1)/2 +1;
Ib(b:b+is-1,b:b+js-1)= image_resize;
image_resize=zeros(size(image_resize));
for i = 1:is
    for j = 1:js
        image_resize(i,j)=sum(sum(Ib(i:i+N-1,j:j+N-1).*h,'omitnan'));
    end
end
axes(handles.axes2);
imshow(image_resize);
a2=image_resize;
% --- Executes on button press in Medianpushbutton.
function Medianpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Medianpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img = a2;
modimg= zeros(size(img)+2);
Outimg = zeros(size(img));
%COPY THE ORIGINAL IMAGE MATRIX TO THE PADDED MATRIX
        for x=1:size(img,1)
            for y=1:size(img,2)
                modimg(x+1,y+1)=img(x,y);
            end
        end    
        
        for i= 1:size(modimg,1)-2
            for j=1:size(modimg,2)-2
                window=zeros(9,1);
                inc=1;
                for x=1:3
                    for y=1:3
                        window(inc)=modimg(i+x-1,j+y-1);
                        inc=inc+1;
                    end
                end 
                med=sort(window);
                Outimg(i,j)=med(5);
            end
        end
%CONVERT THE OUTPUT MATRIX TO 0-255 RANGE IMAGE TYPE
Outimg=uint8(Outimg);

axes(handles.axes2);
imshow(Outimg);
a2=Outimg;
% --- Executes on button press in Lappushbutton.
function Lappushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Lappushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img=a2;
I=zeros(size(img));

%Filter Masks
F1=[1 1 1;1 -8 1; 1 1 1];
%Pad array with zeros
img=padarray(img,[1,1]);
img=double(img);

for i=1:size(img,1)-2
    for j=1:size(img,2)-2
        I(i,j)=sum(sum(F1.*img(i:i+2,j:j+2)));
    end
end
I=uint8(I);
axes(handles.axes2);
imshow(I);
a2=I;
% --- Executes on button press in Sharppushbutton.
function Sharppushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Sharppushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
A=a2;
I1=A;
I=zeros(size(A));
%Filter Masks
F1=[0 1 0;1 -4 1; 0 1 0];
F2=[1 1 1;1 -8 1; 1 1 1];
%Pad array with zeros
A=padarray(A,[1,1]);
A=double(A);
%Implementation of the equation in 2nd Dervitive of laplacian in both
%directions 
for i=1:size(A,1)-2
    for j=1:size(A,2)-2
        I(i,j)=sum(sum(F1.*A(i:i+2,j:j+2)));
    end
end
I=uint8(I);
%Sharpenend Image
B=I1-I;

axes(handles.axes2);
imshow(B);
a2=B;
% --- Executes on button press in Convpushbutton.
function Convpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Convpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
In= a2;
In=double(In);

mask1=[1, 0, -1;1, 0, -1;1, 0, -1]; 
mask2=[1, 1, 1;0, 0, 0;-1, -1, -1]; 
mask3=[0, -1, -1;1, 0, -1;1, 1, 0]; 
mask4=[1, 1, 0;1, 0, -1;0, -1, -1];

mask1=flipud(mask1); 
mask1=fliplr(mask1); 
mask2=flipud(mask2); 
mask2=fliplr(mask2); 
mask3=flipud(mask3); 
mask3=fliplr(mask3); 
mask4=flipud(mask4); 
mask4=fliplr(mask4); 

for i=2:size(In, 1)-1
    for j=2:size(In, 2)-2
        neighbour_matrix1=mask1.*In(i-1:i+1, j-1:j+1); 
        avg_value1=sum(neighbour_matrix1(:)); 
  
        neighbour_matrix2=mask2.*In(i-1:i+1, j-1:j+1); 
        avg_value2=sum(neighbour_matrix2(:)); 
  
        neighbour_matrix3=mask3.*In(i-1:i+1, j-1:j+1); 
        avg_value3=sum(neighbour_matrix3(:)); 
  
        neighbour_matrix4=mask4.*In(i-1:i+1, j-1:j+1); 
        avg_value4=sum(neighbour_matrix4(:)); 
  
        %using max function for detection of final edges 
        In2(i, j)=max([avg_value1, avg_value2, avg_value3, avg_value4]); 
  
    end 
end 
axes(handles.axes2);
imshow(In2);
a2=In2;
% --- Executes on button press in Invpushbutton.
function Invpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Invpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2 ;
in=a2;
for r=1:size(a2,1)
    for c=1:size(a2,2)
        out(r,c,:)=255-a2(r,c,:);
    end
end
axes(handles.axes2);
imshow(out);
a2=out;
% --- Executes on button press in Picpushbutton.
function Picpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Picpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
global a2;
global a3;

a =uigetfile();
filename=a;
setappdata(0,'filename',filename);
a=imread(a);
a3=a;
a = rgb2gray(a);
a2=a;
axes(handles.axes1);
imshow(a);
axes(handles.axes2);
imshow(a2);
axes(handles.axes3);
imshow(a3);
% --- Executes on button press in ChPicpushbutton.
function ChPicpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ChPicpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
global a2;
a =uigetfile();
filename=a;
setappdata(0,'filename',filename);
a=imread(a);

a2 =uigetfile();
filename2=a2;
setappdata(0,'filename2',filename2);
a2=imread(a2);

a = rgb2gray(a);
a2=rgb2gray(a2);
axes(handles.axes1);
imshow(a);
axes(handles.axes2);
imshow(a2);
axes(handles.axes3);
imshow(a2);
% --- Executes on button press in Diffpushbutton.
function Diffpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Diffpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
global a2;

img1 = a;
img2 = a2;
for r=1:size(a2,1)
    for c=1:size(a2,2)
        d(r,c,:) = a2(r,c,:) - a(r,c,:);
    end
end

axes(handles.axes2);
imshow(d);
axes(handles.axes3);
imshow(a2);
a2=d;
% --- Executes on button press in Prewittpushbutton.
function Prewittpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Prewittpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img =a2;
img = double(img);
Fimg = zeros(size(img));

% Prewitt Operator Mask 
Mx = [-1 0 1; -1 0 1; -1 0 1]; 
My = [-1 -1 -1; 0 0 0; 1 1 1]; 

for i = 1:size(img, 1) - 2 
    for j = 1:size(img, 2) - 2 
  
        % Gradient approximations 
        Gx = sum(sum(Mx.*img(i:i+2, j:j+2))); 
        Gy = sum(sum(My.*img(i:i+2, j:j+2))); 
                 
        % Calculate magnitude of vector 
        Fimg(i+1, j+1) = sqrt(Gx.^2 + Gy.^2); 
         
    end
end
% Define a threshold value 
thresholdValue = 100; % varies between [0 255] 
Outimg = max(Fimg, thresholdValue); 
Outimg(Outimg == round(thresholdValue)) = 0; 
Outimg = im2bw(Outimg);

axes(handles.axes2);
imshow(Outimg);
a2=Outimg;
% --- Executes on button press in Sobelpushbutton.
function Sobelpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Sobelpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img = a2;
img = double(img);

for i=1:size(img,1)-2
    for j=1:size(img,2)-2
        %Sobel mask for x-direction:
        SobelX=((2*img(i+2,j+1)+img(i+2,j)+img(i+2,j+2))-(2*img(i,j+1)+img(i,j)+img(i,j+2)));
        %Sobel mask for y-direction:
        SobalY=((2*img(i+1,j+2)+img(i,j+2)+img(i+2,j+2))-(2*img(i+1,j)+img(i,j)+img(i+2,j)));
     
        %The gradient of the image
        Outimg(i,j)=sqrt(SobelX.^2+SobalY.^2);
     
    end
end

axes(handles.axes2);
imshow(Outimg);
a2=Outimg;
% --- Executes on button press in SPNpushbutton.
function SPNpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SPNpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img = a2;
img = imnoise(img , 'salt & pepper',0.02);
axes(handles.axes2);
imshow(img);
a2=img;
% --- Executes on button press in Rconpushbutton.
function Rconpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Rconpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a2;
img = a2;
img2 = img * 0.3 ;
axes(handles.axes2);
imshow(img2);
a2=img2;
% --- Executes on button press in wienerpushbutton.
function wienerpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to wienerpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global a;
img = a;
J = imnoise(img,'gaussian',0,0.025);
K = wiener2(J,[5 5]);
axes(handles.axes2);
imshow(K);
a2=K;
