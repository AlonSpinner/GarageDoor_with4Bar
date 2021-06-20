function varargout = Door_Lift(varargin)
% Door_Lift MATLAB code for Door_Lift.fig
%      Door_Lift, by itself, creates a new Door_Lift or raises the existing
%      singleton*.
%
%      H = Door_Lift returns the handle to a new Door_Lift or the handle to
%      the existing singleton*.
%
%      Door_Lift('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Door_Lift.M with the given input arguments.
%
%      Door_Lift('Property','Value',...) creates a new Door_Lift or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Door_Lift_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Door_Lift_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Door_Lift

% Last Modified by GUIDE v2.5 21-Feb-2019 17:36:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Door_Lift_OpeningFcn, ...
                   'gui_OutputFcn',  @Door_Lift_OutputFcn, ...
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
% --- Executes just before Door_Lift is made visible.
function Door_Lift_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Door_Lift (see VARARGIN)

% Choose default command line output for Door_Lift
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes Door_Lift wait for user response (see UIRESUME)
% uiwait(handles.Fig);
% --- Outputs from this function are returned to the command line.

%build room
ResetRmPush_Callback(hObject, eventdata, handles)
%calculate data points
CalcPush_Callback(hObject, eventdata, handles)
function varargout = Door_Lift_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

%Define varagout as handles
varargout{1} = handles.Fig; 
%% Callbacks
function CalcPush_Callback(hObject, eventdata, handles)
% hObject    handle to CalcPush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Obtain Number of Points for drawings
N=str2num(handles.NPntsEdit.String);
if isempty(N)||N<=0||(mod(N,1)~=0), errordlg('Please input a correct number of points','Four Fingers'); return, end

%Obtain room parameters from RmTable
RmParData=cellfun(@str2num, handles.RmParTable.Data,'un',0);
if any(isempty(RmParData)), errordlg('Some room parameterss werent set properlly','Four Fingers'); return; end
Hroom=RmParData{1}/1e3; Lroom=RmParData{2}/1e3; Hdoor=RmParData{3}/1e3; Hshelf=RmParData{4}/1e3; Hmotor=RmParData{6}/1e3; Mdoor=RmParData{7}; %m

%Obtain data from VarTable
VarData=cellfun(@str2num, handles.VarTable.Data,'un',0);
if any(isempty(VarData)), errordlg('Some variables werent set properlly','Four Fingers'); return; end
Acthd=VarData{1}/1e3; Outhd=VarData{2}/1e3; ActHc=VarData{3}/1e3; OutHc=VarData{4}/1e3; MtrAngleRange=deg2rad(VarData{5});
RFMtr=VarData{6}; RFAct=VarData{7}; RunTime=VarData{8};

%Obtain elbow shapes
RelbowPop=handles.RelbowPop; LelbowPop=handles.LelbowPop;
Relbow=RelbowPop.String{RelbowPop.Value};
Lelbow=LelbowPop.String{LelbowPop.Value};

%R4Bar - find Rp0,Rp1,Ra,Rb and Rf
ActArcP1=Lroom+1i*Acthd; %point in room of actuated axis on door when closed
ActArcP2=(Hdoor-Acthd)+1i*Hshelf; %point in room of actuated axis on door when stored
[Rp0,Ra]=TwoPntHSynth(ActArcP1,ActArcP2,ActHc);
OutArcP1=Lroom+1i*Outhd; %point in room of out axis on door when closed
OutArcP2=(Hdoor-Outhd)+1i*Hshelf; %point in room of out axis on door when stored
[Rp1,Rb]=TwoPntHSynth(OutArcP1,OutArcP2,OutHc);
Rf=Acthd-Outhd;
%L4Bar - find Lp0,Lp1,ActAngle,La,Lf,Lb
Lp0=Hmotor*1i; %m
Lp1=Rp0; %n
ActAngleRange=[angle(ActArcP1-Rp0),angle(ActArcP2-Rp0)]; %rad
[La,Lf,Lb]=L4barSynth(MtrAngleRange,ActAngleRange,RFMtr,RFAct,Lp0,Lp1); %m

%------------Draw 
% DrawDoorArc(Ax,Acthd,R4ActHc,Lroom,Hshelf,Hdoor);
% DrawDoorArc(Ax,Outhd,R4OutHc,Lroom,Hshelf,Hdoor);

%Build Time Vector and MtrAngleVec
MtrVel=(MtrAngleRange(2)-MtrAngleRange(1))/RunTime; %Rad/s
if isempty(MtrVel), errordlg('Please input a correct motor velocity','Four Fingers'); return, end
t=BuildTime(MtrAngleRange(1),MtrAngleRange(2),MtrVel,N);
MtrAngle=linspace(MtrAngleRange(1),MtrAngleRange(2),N);

%Obtain direction and applay on angle vec
Drctn=handles.DirectionPop.Value; %1 to store door, 2 to close door
if Drctn==2
    MtrAngle=fliplr(MtrAngle); 
    MtrVel=-MtrVel;
end

%---------initalize values for graphs
Rf_DoorCG=Acthd-Hdoor/2; %m
Rf_TopEdgeDoor=Acthd-Hdoor; %m
Rf_BotEdgeDoor=Acthd; %m
AxisPnts=zeros(7,length(t));
DoorEdgePnts=zeros(2,length(t));
[DoorCG,DoorCGvel,DoorCGacce,DoorOmega,DoorAlpha,...
    PushAngle,ActAngle,DoorAngle,OutAngle,...
    TPwr,TPwrI,TNw2,TNw2I,FNw2,FNw2I]=deal(zeros(size(t))); %row vectors

%Initalize waitbar
WaitH=waitbar(0,'Calculating Data. Please Wait. . .','Name','Four Fingers');

%Do the loop - calculate all the data and plot it as waiting bar
for kk=1:N
    %Calculate 4bars
    %L4Bar
    [LRmPnts,Lalpha,Lbeta,Lgama]=Calc4BarPnts(Lp0,Lp1,La,Lf,Lb,MtrAngle(kk),Lelbow);
    %R4Bar
    ActAngle(kk)=angle(LRmPnts(3)-LRmPnts(4));
    [RRmPnts,Ralpha,Rbeta,Rgama]=Calc4BarPnts(Rp0,Rp1,Ra,Rf,Rb,ActAngle(kk),Relbow);
    
    %Angular velocities and accelerations
    [Ldgama,~,Ld2gama,~]=AngleVnA(Lalpha,MtrVel,0,Lbeta,Lgama,La,Lf,Lb); %L4bar
    Rdalpha=Ldgama; Rd2alpha=Ld2gama;
    [~,Rdbeta,~,Rd2beta]=AngleVnA(Ralpha,Rdalpha,Rd2alpha,Rbeta,Rgama,Ra,Rf,Rb); %R4bar
    %--------------Calculate graph numbers
    %AxisPnts
    AxisPnts(:,kk)=[LRmPnts(1:3),RRmPnts]; %m 7Xkk
    %Door kinematics
    DoorEdgePnts(:,kk)=[fpnt(RRmPnts(2),RRmPnts(3),Rf_BotEdgeDoor),fpnt(RRmPnts(2),RRmPnts(3),Rf_TopEdgeDoor)]; %m
    DoorCG(:,kk)=mean(DoorEdgePnts(:,kk)); %m
    DoorCGvel(kk)=fVel(Ralpha,Rbeta,Rdalpha,Rdbeta,Rp1,Rp0,Ra,Rf_DoorCG); %m/s
    DoorCGacce(kk)=fAcce(Ralpha,Rbeta,Rdalpha,Rd2alpha,Rdbeta,Rd2beta,Rp1,Rp0,Ra,Rf_DoorCG); %m/s^2
    DoorOmega(kk)=Rdbeta; %rad/s
    DoorAlpha(kk)=Rd2beta; %rad/s^2
    %Motor Torque by Power
    TPwr(kk)=PowerTrq(MtrVel,DoorCGvel(kk),DoorCGacce(kk),DoorOmega(kk),DoorAlpha(kk),Hdoor,Mdoor,0); %Nm
    TPwrI(kk)=PowerTrq(MtrVel,DoorCGvel(kk),DoorCGacce(kk),DoorOmega(kk),DoorAlpha(kk),Hdoor,Mdoor,1); %Nm  
    %Reaction Forces and Torque by Newton2
    PushAngle(kk)=angle(LRmPnts(3)-LRmPnts(2));
    DoorAngle(kk)=angle(RRmPnts(3)-RRmPnts(2));
    OutAngle(kk)=angle(RRmPnts(3)-RRmPnts(4));
    [TNw2(kk),FNw2(kk)]=TrqNw2(MtrAngle(kk),PushAngle(kk),ActAngle(kk),OutAngle(kk),DoorAngle(kk),Ra,Rf,La,Lb,Rf_DoorCG,Hdoor,Mdoor,DoorCGacce(kk),DoorAlpha(kk),0);
    [TNw2I(kk),FNw2I(kk)]=TrqNw2(MtrAngle(kk),PushAngle(kk),ActAngle(kk),OutAngle(kk),DoorAngle(kk),Ra,Rf,La,Lb,Rf_DoorCG,Hdoor,Mdoor,DoorCGacce(kk),DoorAlpha(kk),1);
    
    %update waitbar
    waitbar(kk/N,WaitH,'Calculating Data. Please Wait. . .','Name','Four Fingers');
end

%------------store handles and graphing vectors
%CalcPush - store for PlayPush callback
CalcPush=handles.CalcPush;
CalcPush.UserData{1}=t; %s
CalcPush.UserData{2}=AxisPnts; %m
CalcPush.UserData{3}=N; %[-]
CalcPush.UserData{4}=Hdoor; %m
CalcPush.UserData{5}=Acthd; %m
CalcPush.UserData{6}=Outhd; %m

%Door CG kinematics
DoorCGPush=handles.DoorCGPush;
DoorCGPush.UserData{1}=t; %s
DoorCGPush.UserData{2}=real(DoorCGvel); %m/s Xvel
DoorCGPush.UserData{3}=imag(DoorCGvel); %m/s Yvel
DoorCGPush.UserData{4}=abs(DoorCGacce); %m/s^2 Xacce
DoorCGPush.UserData{5}=DoorAlpha; %rad/s^2 Yacce
%Motor Torque
MotorTorquePush=handles.MotorTorquePush;
MotorTorquePush.UserData{1}=t;
MotorTorquePush.UserData{2}=TPwr;
MotorTorquePush.UserData{3}=TNw2;
MotorTorquePush.UserData{4}=TPwrI;
MotorTorquePush.UserData{5}=TNw2I;
%Check Collision
CollisionPush=handles.CollisionPush;
CollisionPush.UserData{1}=t; %[s]
CollisionPush.UserData{2}=DoorEdgePnts; %[m]
CollisionPush.UserData{3}=AxisPnts(2,:); %second axis point %[m]
CollisionPush.UserData{4}=Hroom; %[m]
%Assignin (store to save to workspace)
AssigninPush=handles.AssigninPush;
AssigninPush.UserData{1}=t; %time vec. 1xN  real [s]
AssigninPush.UserData{2}=AxisPnts; %6XN complex [m]
AssigninPush.UserData{3}=DoorEdgePnts; %2XN complex [m]
AssigninPush.UserData{4}=DoorCG; %1XN complex [m]
AssigninPush.UserData{5}=abs(DoorCGvel); %1XN real [m/s]
AssigninPush.UserData{6}=DoorOmega; %1XN real [rad/s]
AssigninPush.UserData{7}=TPwr; %1XN real [Nm]
AssigninPush.UserData{8}=TPwrI; %1XN real [Nm]
AssigninPush.UserData{9}=MtrVel; %scalar real [rad/s]

%Re-action forces
ReactionPush=handles.ReactionPush;
ReactionPush.UserData{1}=t;
ReactionPush.UserData{2}=FNw2;
ReactionPush.UserData{3}=FNw2I;

%turn PlotButtons green - available for plotting
GrassGreen=[0.47,0.67,0.1];
handles.DoorCGPush.BackgroundColor=GrassGreen;
handles.ReactionPush.BackgroundColor=GrassGreen;
handles.MotorTorquePush.BackgroundColor=GrassGreen;
handles.AssigninPush.BackgroundColor=GrassGreen;
handles.CollisionPush.BackgroundColor=GrassGreen;
handles.PlayPush.BackgroundColor=GrassGreen;
handles.MoviePush.BackgroundColor=GrassGreen;

%kill waitbar
delete(WaitH);
function StopPush_Callback(hObject, eventdata, handles)
%1 is stop, 0 is contiue
if hObject.UserData
    hObject.UserData=0;
else
    hObject.UserData=1;
end
function ResetRmPush_Callback(hObject, eventdata, handles)
Ax=handles.Ax;
cla(Ax,'reset');

%turn PlotPush buttons red - won't work
LightRed=[0.9,0.6,0.6];
handles.DoorCGPush.BackgroundColor=LightRed;
handles.ReactionPush.BackgroundColor=LightRed;
handles.MotorTorquePush.BackgroundColor=LightRed;
handles.AssigninPush.BackgroundColor=LightRed;
handles.CollisionPush.BackgroundColor=LightRed;
handles.PlayPush.BackgroundColor=LightRed;
handles.MoviePush.BackgroundColor=LightRed;

%Obtain room parameters from RmTable
RmParData=cellfun(@str2num, handles.RmParTable.Data,'un',0);
if any(isempty(RmParData)), errordlg('Some room parameterss werent set properlly','Four Fingers'); return; end
Hroom=RmParData{1}; Lroom=RmParData{2}; Hdoor=RmParData{3}; Hshelf=RmParData{4}; Lshelf=RmParData{5}; Hmotor=RmParData{6};

%Obtain data from VarTable
VarData=cellfun(@str2num, handles.VarTable.Data,'un',0);
if any(isempty(VarData)), errordlg('Some variables werent set properlly','Four Fingers'); return; end
Acthd=VarData{1}; Outhd=VarData{2}; ActHc=VarData{3}; OutHc=VarData{4}; MtrAngle=deg2rad(VarData{5});
RFMtr=VarData{6}; RFAct=VarData{7};

%Obtain elbow shapes
RelbowPop=handles.RelbowPop; LelbowPop=handles.LelbowPop;
Relbow=RelbowPop.String{RelbowPop.Value};
Lelbow=LelbowPop.String{LelbowPop.Value};

%---------do some initial drawing
Ax=handles.Ax;
DrawRoom(Ax,Hroom,Lroom,Hdoor,Hshelf,Lshelf,Hmotor)
door_handle=DrawDoor(Ax,Hdoor,Lroom+1i*Acthd,Lroom+1i*Outhd,Acthd,Outhd);

%------------calculate R4bar
%find Lp0,Lp1,La,Lb and Lf
ActArcP1=Lroom+1i*Acthd; %point in room of actuated axis on door when closed
ActArcP2=(Hdoor-Acthd)+1i*Hshelf; %point in room of actuated axis on door when stored
[Rp0,Ra]=TwoPntHSynth(ActArcP1,ActArcP2,ActHc);
OutArcP1=Lroom+1i*Outhd; %point in room of out axis on door when closed
OutArcP2=(Hdoor-Outhd)+1i*Hshelf; %point in room of out axis on door when stored
[Rp1,Rb]=TwoPntHSynth(OutArcP1,OutArcP2,OutHc);
Rf=Acthd-Outhd;
%Calc points
R4bar_RmPnts=Calc4BarPnts(Rp0,Rp1,Ra,Rf,Rb,angle(ActArcP1-Rp0),Relbow);

%------------calculate L4bar
%find La,Lf,Lb
Lp0=Hmotor*1i;
Lp1=Rp0;
ActAngle=[angle(ActArcP1-Rp0),angle(ActArcP2-Rp0)];
[La,Lf,Lb]=L4barSynth(MtrAngle,ActAngle,RFMtr,RFAct,Lp0,Lp1);
%Calc points
L4bar_RmPnts=Calc4BarPnts(Lp0,Lp1,La,Lf,Lb,MtrAngle(1),Lelbow);

%Draw 4bar mechanisms
L4bar_handle=DrawFourBar(Ax,L4bar_RmPnts);
R4bar_handle=DrawFourBar(Ax,R4bar_RmPnts);

%Draw Arcs
% NicePurp=[0.7,0,0.7];
% DrawDoorArc(Ax,Acthd,ActHc,Lroom,Hshelf,Hdoor,'b');
% DrawDoorArc(Ax,Outhd,OutHc,Lroom,Hshelf,Hdoor,NicePurp);

%store drawing handles in PlayPush user data
handles.PlayPush.UserData=[door_handle,R4bar_handle,L4bar_handle];

%stores links length BarLinksPush user data
Lg=abs(Lp1-Lp0); Rg=abs(Rp1-Rp0);
handles.BarLinksPush.UserData=[La,Lf,Lb,Lg,Ra,Rf,Rb,Rg];
function PlayPush_Callback(hObject, eventdata, handles)
LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%delete old drawings from user data
delete(handles.PlayPush.UserData);

%obtain pausetime, Ax, and calculated data (change m->mm in process)
Ax=handles.Ax;
PauseTime=str2num(handles.PauseTimeEdit.String);
if isempty(PauseTime)||PauseTime<=0, errordlg('Please input a correct PauseTime','Four Fingers'); return, end
CalcPush=handles.CalcPush;
t=CalcPush.UserData{1}; %s
AxisPnts=CalcPush.UserData{2}*1e3; %mm
N=CalcPush.UserData{3}; %[-]
Hdoor=CalcPush.UserData{4}*1e3; %mm
Acthd=CalcPush.UserData{5}*1e3; %mm
Outhd=CalcPush.UserData{6}*1e3; %mm

%Initalize handles for drawing
[R4bar_handle,door_handle,L4bar_handle,Time_handle]=deal([]);

%Obtain StopPush. LightRed means go. DarkRed means stop
StopPush=handles.StopPush;

%split data and turn to mm
LPnts=AxisPnts(1:4,:); %mm
RPnts=AxisPnts(4:7,:); %mm

%Do the plotting - simulation
for kk=1:N
    if StopPush.UserData %if Stop was pressed during run
       StopPush.UserData=0; %unstop
        return
    end 
    delete([R4bar_handle,door_handle,L4bar_handle,Time_handle]);
    
    %Drawings
    L4bar_handle=DrawFourBar(Ax,LPnts(:,kk)); %turn to mm
    door_handle=DrawDoor(Ax,Hdoor,RPnts(2,kk),RPnts(3,kk),Acthd,Outhd); %turn to mm
    R4bar_handle=DrawFourBar(Ax,RPnts(:,kk));%turn to mm
    Time_handle=DrawTimeStamp(Ax,t(kk));
    pause(PauseTime); %as load time
end

%Drawings in CalcPush
handles.PlayPush.UserData=[door_handle,R4bar_handle,L4bar_handle,Time_handle];
function MoviePush_Callback(hObject, eventdata, handles)
LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%Obtain room parameters from RmTable
RmParData=cellfun(@str2num, handles.RmParTable.Data,'un',0);
if any(isempty(RmParData)), errordlg('Some room parameterss werent set properlly','Four Fingers'); return; end
Hroom=RmParData{1}; Lroom=RmParData{2}; Hdoor=RmParData{3}; Hshelf=RmParData{4}; Lshelf=RmParData{5}; Hmotor=RmParData{6}; %sizes in mm

%Create new invisible ax and draw room on them
Fig=figure('visible','off');
Ax=axes(Fig);
DrawRoom(Ax,Hroom,Lroom,Hdoor,Hshelf,Lshelf,Hmotor)

%obtain movie name, and calculated data (change m->mm in process)
MovieName=handles.MovieEdit.String;
CalcPush=handles.CalcPush;
t=CalcPush.UserData{1}; %s
AxisPnts=CalcPush.UserData{2}*1e3; %mm
N=CalcPush.UserData{3}; %[-]
Acthd=CalcPush.UserData{5}*1e3; %mm
Outhd=CalcPush.UserData{6}*1e3; %mm

%Initalize handles for drawing
[R4bar_handle,door_handle,L4bar_handle,Time_handle]=deal([]);

%Obtain StopPush. LightRed means go. DarkRed means stop
StopPush=handles.StopPush;

%split data and turn to mm
LPnts=AxisPnts(1:4,:); %mm
RPnts=AxisPnts(4:7,:); %mm

%Initailize video file
vidWriter = VideoWriter([MovieName,'.mp4'],'MPEG-4');
vidWriter.FrameRate=round(N/t(end));
open(vidWriter);

%Initalize waitbar
WaitH=waitbar(0,'Rendering Movie. Please Wait. . .','Name','Four Fingers');

%Do the plotting - simulation
for kk=1:N
    if StopPush.UserData %if Stop was pressed during run
       StopPush.UserData=0; %unstop
        return
    end 
    delete([R4bar_handle,door_handle,L4bar_handle,Time_handle]);
    
    %Drawings
    L4bar_handle=DrawFourBar(Ax,LPnts(:,kk)); %turn to mm
    door_handle=DrawDoor(Ax,Hdoor,RPnts(2,kk),RPnts(3,kk),Acthd,Outhd); %turn to mm
    R4bar_handle=DrawFourBar(Ax,RPnts(:,kk));%turn to mm
    Time_handle=DrawTimeStamp(Ax,t(kk));
    writeVideo(vidWriter,getframe(Ax));
    
    %update waitbar
    waitbar(kk/N,WaitH,'Rendering Data. Please Wait. . .','Name','Four Fingers');
end

%Drawings in CalcPush
handles.PlayPush.UserData=[door_handle,R4bar_handle,L4bar_handle,Time_handle];

%Close VidWriter
close(vidWriter);

%kill WaitBar
close(WaitH);

%inform user
helpdlg('Movie saved to matlab working folder','Four Fingers');

%plotting callbacks      
function DoorCGPush_Callback(hObject, eventdata, handles)
%returns graphs for door CG motion:
%X Velocity(t), Y Velocity(t), Abs Acceleration(t), Radial Accleeration(t)

LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%------------Obtain vectors of interest
t=hObject.UserData{1};
Xvel=hObject.UserData{2};
Yvel=hObject.UserData{3};
Absacce=hObject.UserData{4};
Alpha=hObject.UserData{5};

%Plot it up
LightBrown=[0.93,0.9,0.6];
Colors=lines(4);

Fig1=figure('name','Door CG kinematics -VelX, VelY ','color',LightBrown);
Fig2=figure('name','Door CG kinematics - Abs Acceleration','color',LightBrown);
Fig3=figure('name','Door CG kinematics - Radial Acceleration','color',LightBrown);
Ax1=axes('parent',Fig1); %X direction Velocity
Ax2=axes('parent',Fig2); %Y direction Velocity
Ax3=axes('parent',Fig3); %X direction Acceleration
grid(Ax1,'on'); grid(Ax2,'on'); grid(Ax3,'on');
hold(Ax1,'on'); hold(Ax2,'on'); hold(Ax3,'on');
plot(Ax1,t,Xvel,'color',Colors(1,:),'linew',2);
plot(Ax1,t,Yvel,'color',Colors(2,:),'linew',2);
plot(Ax2,t,Absacce,'color',Colors(3,:),'linew',2);
plot(Ax3,t,Alpha,'color',Colors(4,:),'linew',2);

xlabel(Ax1,'[s]'); ylabel(Ax1,'[m/s]');
xlabel(Ax2,'[s]'); ylabel(Ax2,'[m/s^2]');
xlabel(Ax3,'[s]'); ylabel(Ax3,'[rad/s^2]');
title(Ax1,'X,Y velocity vs Time ');
title(Ax2,'Absolute Acceleration vs Time ');
title(Ax3,'Radial Acceleration vs Time ');
legend(Ax1,'X Velocity','Y Velocity');
function MotorTorquePush_Callback(hObject, eventdata, handles)
%returns Motor Torque graphing

LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%------------Obtain vectors of interest
t=hObject.UserData{1};
TPwr=hObject.UserData{2};
TNw2=hObject.UserData{3};
TPwrI=hObject.UserData{4};
TNw2I=hObject.UserData{5};

%Plot it up
LightBrown=[0.93,0.9,0.6];
Colors=lines(2);

Fig=figure('name','Motor Torque ','color',LightBrown);
Ax=axes('parent',Fig); %without interia
grid(Ax,'on'); hold(Ax,'on');
plot(Ax,t,TPwr,'linew',2);
plot(Ax,t,TNw2,'linew',2);
plot(Ax,t,TPwrI,'linew',2,'color',Colors(1,:),'lines','--');
plot(Ax,t,TNw2I,'linew',2,'color',Colors(2,:),'lines','--');
xlabel(Ax,'[s]'); ylabel(Ax,'Nm');
title(Ax,'Motor Torque vs Time');
legend(Ax,'Power','Nw2','Power /w Inertia','Nw2 /w Inertia','location','best');
function ReactionPush_Callback(hObject, eventdata, handles)
%plots reaction forces in motor axis

LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%------------Obtain vectors of interest
t=hObject.UserData{1}; %s
F=hObject.UserData{2}; %N without inertia
FI=hObject.UserData{3}; %N with inertia

%Plot it up
LightBrown=[0.93,0.9,0.6];
Colors=lines(3);

Fig=figure('name','Reaction Forces in Motor Axis ','color',LightBrown);
Ax1=axes('parent',Fig);
grid(Ax1,'on');
hold(Ax1,'on');
plot(Ax1,t,real(F),'linew',2); %X axis
plot(Ax1,t,imag(F),'linew',2); %Y axis
plot(Ax1,t,abs(FI),'linew',2); %abs
plot(Ax1,t,real(FI),'linew',2,'color',Colors(1,:),'lines','--'); %X axis
plot(Ax1,t,imag(FI),'linew',2,'color',Colors(2,:),'lines','--'); %Y axis
plot(Ax1,t,abs(FI),'linew',2,'color',Colors(3,:),'lines','--'); %abs
xlabel(Ax1,'[s]'); ylabel(Ax1,'[N]');
title(Ax1,'Reaction Forces in Motor Axis vs Time');
legend(Ax1,'X','Y','abs','X /w Inertia','Y /w Inertia','abs /w Inertia');
function CollisionPush_Callback(hObject, eventdata, handles)
LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please calculate data prior','Four Fingers');
    return
end

%------------Obtain vectors of interest
t=hObject.UserData{1};
TopEdge=hObject.UserData{2}(1,:); %m
BotEdge=hObject.UserData{2}(2,:); %m
PushAxis=hObject.UserData{3}; %m
Hroom=hObject.UserData{4}; %m

%Plot it up
LightBrown=[0.93,0.9,0.6];

Fig=figure('name','Door Collision','color',LightBrown);
Ax=axes('parent',Fig); %without interia
grid(Ax,'on'); hold(Ax,'on');
plot(Ax,t,imag(TopEdge),'linew',2);
plot(Ax,t,imag(BotEdge),'linew',2);
plot(Ax,t,imag(PushAxis),'linew',2);
plot(Ax,Ax.XLim,[Hroom,Hroom],'linew',2,'color','k');
plot(Ax,Ax.XLim,[0,0],'linew',2,'color','k');
xlabel(Ax,'[s]'); ylabel(Ax,'Height [m]');
title(Ax,'Point Height vs Time');
legend(Ax,'Top Door Edge','Bot Door Edge','Push Axis','location','best');

%info callbacks
function BarLinksPush_Callback(hObject, eventdata, handles)
%show user the links length of the left four bar mechanism
%Obtain room parameters from RmTable

LinksLength=hObject.UserData;

%show user the links length of the left four bar mechanism
helpdlg(sprintf(['Left Four Bar Links Length [mm]\n\n',...
    'a = %g\nf = %g\nb = %g\ng = %g\n\n',...
    'Right Four Bar Links Length [mm]\n\n',...
    'a = %g\nf = %g\nb = %g\ng = %g'],LinksLength),'Four Fingers');
function AssigninPush_Callback(hObject, eventdata, handles)
LightRed=[0.9,0.6,0.6];
if norm(hObject.BackgroundColor-LightRed)<eps
    errordlg('Please play the mechanism before continuing','Four Fingers');
    return
end

%prints vectors to base workspace
t=hObject.UserData{1}; %time vec. 1xN  real [s]
AxisPnts=hObject.UserData{2}; %7XN complex [m]
DoorEdgePnts=hObject.UserData{3}; %2XN complex [m]
DoorCG=hObject.UserData{4}; %1XN complex [m]
absDoorCGvel=hObject.UserData{5}; %1XN complex [m/s]
DoorOmega=hObject.UserData{6}; %1XN real [rad/s]
PwrMtrTrq=hObject.UserData{7}; %1XN real [Nm]
PwrIntMtrTrq=hObject.UserData{8}; %1XN real [Nm]
MtrVel=hObject.UserData{9}; %scalar real [rad/s]

assignin('base','Tvec',t); %time vec. 1xN  real [s]
assignin('base','r_i',AxisPnts); %7XN complex [m]
assignin('base','p_i',DoorEdgePnts); %2XN complex [m]
assignin('base','rc',DoorCG); %1XN complex [m]
assignin('base','vc',absDoorCGvel); %1XN complex [m/s]
assignin('base','w',DoorOmega);  %1XN real [rad/s]
assignin('base','tau0',PwrMtrTrq); %1XN real [Nm]
assignin('base','tauI',PwrIntMtrTrq);  %1XN real [Nm]
assignin('base','omega',MtrVel); %scalar real [rad/s]
%% Functions 
function nothing %just to make things with order. allows for Fourbar title

%Fourbar computation
function delta=CalcDelta(g,a,f,b,alpha)
%Input: 
%g - ground
%a - actuator link
%f - middle link
%b - output angle link
%alpha - angle of actuator link (entrance), radians
%p0 - point where a meets g
%p1 - point where b meets g

%Output:
%delta: rad. returns real(delta)

r0=g; %r ~ direction and size vector of link
r1=a*exp(1i*alpha);
Px=real(r1-r0);
Py=imag(r1-r0);
delta=acos((b^2+f^2-Px^2-Py^2)/(2*f*b)); %rad
delta=real(delta); %to fix numerical issues
function [gama,beta]=CalcGamaBeta(g,a,f,b,alpha,delta)
%Used in fcn Calc4BarPnts

%g - ground
%a - actuator link
%f - middle link
%b - output angle link
%alpha - angle of actuator link (entrance), radians
%delta-rad

r0=g; %r ~ direction and size vector of link
r1=a*exp(1i*alpha);
Px=real(r1-r0);
Py=imag(r1-r0);

A=[b*cos(delta)-f,-b*sin(delta);
    b*sin(delta),b*cos(delta)-f];
b=[Px;Py];
x=A\b; %[cos(beta);sin(beta)]W

beta=real(atan2(x(2),x(1))); %to fix numerical issues
gama=beta+delta;
gama=real(gama);
function [RmPnts,alpha,beta,gama]=Calc4BarPnts(p0,p1,a,f,b,theta,elbow)
%g - ground
%a - actuator link
%f - middle link
%b - output angle link
%theta - angle (rad) of actuator link corresponding to the horrizon
%alpha - angle (rad) of actuator link corresponding to mechanism ground
%p0 - point where a meets g
%p1 - point where b meets g
%elbow - elbow up or down of solution. string

switch elbow
    case 'Up'
        deltasign=1;
    case 'Down'
        deltasign=-1;
end

%Output:
%points representing axis of links - [p0,aXf,fXb,p1]
alpha=+theta-angle(p1-p0);
g=abs(p1-p0); %ground link length
delta=CalcDelta(g,a,f,b,alpha);
[gama,beta]=CalcGamaBeta(g,a,f,b,alpha,delta*deltasign);
r0=g; %r ~ direction and size vector of link
r1=a*exp(1i*alpha);
r3=b*exp(1i*gama);
HomPnts=[0,r1,r0+r3,r0]; %Points in link system
RmPnts=HomPnts*exp(1i*angle(p1-p0))+p0; %points in room system
function [dgama,dbeta,d2gama,d2beta]=AngleVnA(alpha,dalpha,d2alpha,beta,gama,a,f,b) %By Yam Ben Natan copywrite
% calculate angular velocities and accelerations in 4bar local system (link
% g is parralel to horrizon)

%alpha - angle of "a"X"g" - input angle 
%dalpha - radial velocity (rad/s)
%d2alpha - radial acceleration (rad/s^2)
%beta - angle of "a"X"f"
%g - ground
%a - actuator link
%f - middle link
%b - output angle link

A=a/sin(gama-beta)*[sin(alpha-gama)/f; 
                        sin(alpha-beta)/b];
dvec = dalpha*A;
dbeta = dvec(1);   dgama = dvec(2); 

B = a/(sin(gama-beta))^2.*...
[(cos(alpha-gama)*(dalpha-dgama)*sin(gama-beta)-sin(alpha-gama)*cos(gama-beta)*(dgama-dbeta))/f;
 (cos(alpha-beta)*(dalpha-dbeta)*sin(gama-beta)-sin(alpha-beta)*cos(gama-beta)*(dgama-dbeta))/b];

d2vec = A.*d2alpha+B.*dalpha; 
d2beta = d2vec(1);   d2gama = d2vec(2);

%simulation vector
function t=BuildTime(MtrAngle_0,MtrAngle_f,MtrVel,N)
%MtrAngles are in rad
%MtrSpeed is in rad/s
%N - number of partitions - length(t)

%returns t - row vector of time
t=(linspace(0,(MtrAngle_f-MtrAngle_0)/MtrVel,N));
function fR=fpnt(aXf,fXb,R) %Door Edges and CG
%aXf,fXb - points of axies in room ordinates
%p1,p0 grounding points. p0 is "a"X"g"
%R - distance of point in the direction of f from atip.
%for doorCG: R=ActHd-Hdoor/2
fR=aXf+R*exp(1i*angle((fXb-aXf)));
function fRvel=fVel(alpha,beta,dalpha,dbeta,p1,p0,a,R) %Door CG Velocity
%dalpha - alpha dot. actuating link radial velocity (rad/s)
%dbeta - beta dot. radial speed of link "f" (rad/s)
%p1,p0 grounding points. p0 is "a"X"g"
%R - distance of point in the direction of f from atip.
%for doorCG: R=ActHd-Hdoor/2

%fRvel - velocity of point R along linkage f in room (complex number)
psi=angle(p1-p0); %angle of ground relative to room
Vatip_O=dalpha*a*exp(1i*(alpha+psi+pi/2)); %velocity of point atip in relation to world O
VfR_atip=dbeta*R*exp(1i*(beta+psi+pi/2)); %velocity of point fR in relation to Vatip
fRvel=Vatip_O+VfR_atip; 
function fRacc=fAcce(alpha,beta,dalpha,d2alpha,dbeta,d2beta,p1,p0,a,R) %Door CG Acceleration
%dalpha - alpha dot. actuating link radial velocity (rad/s)
%dbeta - beta dot. radial speed of link "f" (rad/s)
%d2alpha - alpha dotayim. actuating link radial acceleration (rad/s^2)
%d2beta - beta dotayim. radial acceleration of link "f" (rad/s^2)
%p1,p0 grounding points. p0 is "a"X"g"
%R - distance of point in the direction of f from atip.
%for doorCG: R=ActHd-Hdoor/2

%fRvel - velocity of point R along linkage f in room (complex number)

%Theory: for a point with constant distance from its rotational axis
%a_r=omega^2*r in -r^ direction
%a_theta=alpha*r in theta^ direction

psi=angle(p1-p0); %angle of ground relative to room
Aatip_O=dalpha^2*a*exp(1i*(alpha+psi+pi))+d2alpha*a*exp(1i*(alpha+psi+pi/2)); %velocity of point atip in relation to world O
AfR_atip=dbeta^2*R*exp(1i*(beta+psi+pi))+d2beta*R*exp(1i*(beta+psi+pi/2)); %velocity of point fR in relation to Vatip
fRacc=Aatip_O+AfR_atip;
function T=PowerTrq(MtrVel,DoorCGVel,DoorCGAcce,DoorOmega,DoorAlpha,Hdoor,Mdoor,IFactor)
%M - kg
%MtrVel - scalar rad/s
%DoorCGVel - m/s complex number. real for X and imag for Y
%DoorCGAcce - m/s same as DoorCGVel
%Hdoor - m
%DoorOmega - rad/s
%DoorAlpha rad/s^2
%IFactor - boolean 0 or 1. 1 if we want to calcualte with inertia

%T*MtrVel=d/dt(M*g*Y+0.5I*omega^2+0.5*M*Vcg^2)=M*g*VCg_y+I*omega*alpha+M*dot(Vcg,Acg)

g=9.81; %m/s^2
I=(1/12)*Mdoor*Hdoor^2; %kg*m^2
Vy=imag(DoorCGVel);
T=(Mdoor*g*Vy+IFactor*(I*DoorOmega*DoorAlpha+Mdoor*DotCmplx(DoorCGVel,DoorCGAcce)))/MtrVel;
function [T,F]=TrqNw2(Mtr,Push,Act,Out,Door,Ra,Rf,La,Lb,C,Hdoor,M,a,alpha,IFactor)
% right Mechanism Parameters
% Ra, Rf, Rb - m
%Act,Door,Out,Act - angles in rad

% left Mechanism Parameters
% La, Lf, Lb - m
%Mtr, Psh, Act - angles in rad

%C - distance from RaXRf to door CG
% M- Mass of door
% Hdoor - m
%a - m/s^2
%alpha - rad/s^2

%IFactor 1 if interia is to be accounted for. 0 otherwise

g=9.81; %m/s^2
I=(1/12)*M*Hdoor^2; %kg*m^2
ax=real(a); %m/s^2
ay=imag(a); %m/s^2


A=[
    1  -1  0  0  0  0  0;
    0  0  0  1  -1  0  0;
    -Rf*sin(Door)  0  0  Rf*cos(Door)  0  0  0;
    tan(Out)  0  0  -1  0  0  0;
    0  1  1  0  0  0  cos(Push);
    0  0  0  0  1  1  sin(Push);
    0  -sin(Act)*Ra  0  0  cos(Act)*Ra  0  Lb*sin(Push-Act);];
B=[IFactor*M*ax; IFactor*M*ay+M*g; IFactor*I*alpha+M*g*cos(Door)*C; 0; 0; 0; 0;];
    
X=linsolve(A,B);

F=X(end);
T=-La*F*sin(Mtr-Push);
F=F*exp(1i*(Push+pi) );

%Drawing functions
function DrawRoom(Ax,Hroom,Lroom,Hdoor,Hshelf,Lshelf,Hmotor)
RoomColor=[0,0,0];
MotorColor=[1,0,1];
RoomWidth=4;

hold(Ax,'on');
axis(Ax,'equal');
Ax.Color=0.999*ones(3,1);

%Plot walls - in mm, are given
P=[Lroom+Hdoor*1i,Lroom+Hroom*1i,Hroom*1i,Hshelf*1i,Lshelf+Hshelf*1i,Hshelf*1i,0,Lroom];
plot(P,'color',RoomColor,'linew',RoomWidth);

%Plot Motor
r=100;
circ=r*exp(1i*linspace(0,2*pi,10));
plot(Hmotor*1i+circ,'color',MotorColor,'linew',RoomWidth);

Ax.XLim(2)=Ax.XLim(2)*1.05;
axis(Ax,'manual'); %fix axes limits
function h=DrawFourBar(Ax,RmPnts)
%Plots FourBar and returns handle
%RmPnts - points representing axis of links - [p0,aXf,fXb,p1]

% Link0=[RmPnts(1),RmPnts(4)]; %g
Link1=[RmPnts(1),RmPnts(2)]; %a
Link2=[RmPnts(2),RmPnts(3)]; %f
Link3=[RmPnts(3),RmPnts(4)]; %b

% plot(Ax,real(Link0),imag(Link0),'k','linewidth',2); %link0 %ground
NicePurp=[0.7,0,0.7];
h1=plot(Ax,real(Link1),imag(Link1),'b','linewidth',2); %link1
h2=plot(Ax,real(Link2),imag(Link2),'r','linewidth',2); %link2
h3=plot(Ax,real(Link3),imag(Link3),'color',NicePurp,'linewidth',2); %link3

phi=linspace(0,2*pi,10);
RBearing=50;
Circ=RBearing*exp(1i*phi);
h4=plot(Ax,RmPnts(1)+Circ,'-k');
h5=plot(Ax,RmPnts(2)+Circ,'-k');
h6=plot(Ax,RmPnts(3)+Circ,'-k');
h7=plot(Ax,RmPnts(4)+Circ,'-k');

h=[h1,h2,h3,h4,h5,h6,h7];
function h=DrawDoor(Ax,Hdoor,Ptop,Pbot,htop,hbot)
%Hdoor - height of door
%Ptop - top axis point in door
%Pbot - bottom axis point in door
%htop - distance of Ptop from buttom edge of door
%hbot - distancce of Pbot from buttom edoge of door
%Etop - top edge of door in room
%Ebot - buttom edge of door in room

%draws to Ax and returns handle

t=(Ptop-Pbot)/abs(Ptop-Pbot);
Etop=Ptop+t*(Hdoor-htop);
Ebot=Pbot-t*hbot;
GrassGreen=[0.1,0.6,0.1];
h=plot(Ax,[Ebot,Etop],'color',GrassGreen,'linewidth',5);
function h=DrawTimeStamp(Ax,t)
h=text(200,200,sprintf('Time: %g [s]',t),'FontSize',12,'Parent',Ax);
function h=DrawDoorArc(Ax,hd,Hc,Lroom,Hshelf,Hdoor,ArcColor)
%Draws Arcs to decide on first fourbar mechanism

%h - distance of axis point along the door from its buttom
%Hc - height of center of arc - to be axis

%draws door arc and return handle

%find two points on arc
P1=Lroom+1i*hd; %on door when closed
P2=(Hdoor-hd)+1i*Hshelf; %on door when stored
%find center of arc
C=TwoPntHSynth(P1,P2,Hc);
%Plot top arc
R=abs(P1-C);
Angle1=angle(P1-C);
Angle2=angle(P2-C);
Circ=C+R*exp(1i*linspace(Angle1,Angle2,30));
h=plot(Ax,Circ,'color',ArcColor,'lines','--','linew',2);

%FourBar two point Synthasis by height from floor
function [C,R]=TwoPntHSynth(P1,P2,Hc)
%Return Point C which has height Hc from the buttom of the room, and has
%the same distance to P1 and P2. 
%a bar with an axis on C will ensure the movement from P1 to P2

%find center of top arc
PM=mean([P1,P2]);
t=(P2-P1)/abs(P2-P1);
n=t*exp(1i*pi/2);
ImpLine1=Pnts2ImpLine(PM+n,+PM-n); %direction n from PtopM
ImpLine2=Pnts2ImpLine(0+1i*Hc,5000+1i*Hc); %line in height Hc
C=LineX(ImpLine1,ImpLine2);
R=abs(P2-C);
function ImpLine=Pnts2ImpLine(p1,p2)
%Ax+By+C=0
%ImpLine=[A,B,C]
%p1,p2 represent points in complex numbers
p1=[real(p1),imag(p1)]; p2=[real(p2),imag(p2)];
ImpLine=[(p2(2)-p1(2)),-(p2(1)-p1(1)),(p1(2)*p2(1)-p2(2)*p1(1))];
function p=LineX(ImpLine1,ImpLine2)
%Ax+By+C=0
%ImpLine=[A,B,C]
%return p - complex number representing the point

%We solve both equations with linear system M[x;y]=N
M=[ImpLine1(1),ImpLine1(2);
    ImpLine2(1),ImpLine2(2)];
N=[-ImpLine1(3);-ImpLine2(3)];

if abs(det(M))<eps %lines are parallel 
    p=[nan,nan];
    return
else %solve the linear algebra
p=M\N; %returns column vec.
p=p(1)+1i*p(2);
end
%L4Bar synth
function [a,f,b]=L4barSynth(MtrAngle,ActAngle,RFMtr,RFAct,p0,p1)
%MtrAngle - 1x2 numeric vector of angles in rad for range
%ActAngle - 1x2 numeric vector of angles in rad for range
%RFMtr, RFAct - range factors are numbers between 0 and 1.
%when MtrAngle did RFMtr*MtrAngle_range of the way, ActAngle did
%RFAct*ActAngle_range

%change to mechanism system:
alpha=MtrAngle-angle(p1-p0);
gama=ActAngle-angle(p1-p0);

alphaRange=alpha(2)-alpha(1);
gamaRange=gama(2)-gama(1);

alpha3=[alpha(1),alphaRange*RFMtr+alpha(1),alpha(2)];
gama3=[gama(1),gamaRange*RFAct+gama(1),gama(2)];

%Define matrices to solve synthasis
T = [
    cos(gama3(1)),cos(alpha3(1)),1;
    cos(gama3(2)),cos(alpha3(2)),1;
    cos(gama3(3)),cos(alpha3(3)),1;
    ];
B = [
    cos(gama3(1)-alpha3(1));
    cos(gama3(2)-alpha3(2));
    cos(gama3(3)-alpha3(3));
    ];
K=inv(T)*B; %algebric functions connecting ratios between links length

g=abs(p1-p0); %ground link length
b=-g/K(2);
a=g/K(1);
f=sqrt(-K(3)*2*a*b+g^2+a^2+b^2);

%Extra Functions
function DotProduct=DotCmplx(z1,z2)
DotProduct=real(z1)*real(z2)+imag(z1)*imag(z2);

%Useful gui functions
function VarTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VarTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%Set VarTable propeties to fix guide issues
hObject.Data=cell(8,1);
hObject.ColumnName=hObject.ColumnName(1);
hObject.ColumnFormat={[]};
hObject.ColumnWidth={75};
hObject.Position(4)=12;

%insert default raw data is nothing is inserted
VarTableData=cellfun(@num2str, {500; %Actuated Door Axis H
               200; %Outlink Door Axis Height
               0; %R4bar Actuator Axis Height
               300; %R4bar Outlink Axis Height
                [-22,90]; %Motor angle range (deg)
                0.39; %Motor angle range factor
                0.25%actuator range factor
                20},...%run time
                'un',0);

hObject.Data=VarTableData;
function RmParTable_CreateFcn(hObject, eventdata, handles)
%Set RmParTable propeties to fix guide issues
hObject.Data=cell(7,1);
hObject.ColumnName=hObject.ColumnName(1);
hObject.ColumnFormat={[]};
hObject.ColumnWidth={75};
hObject.Position(4)=10.5;

%insert default raw data is nothing is inserted
RmParTableData=cellfun(@num2str, {3000; %Room Height
               5000; %Room Length
               2000; %Door Height
               2300; %Shelf Height
                1500; %Shelf Length
                1500; %Motor Height
                40},'un',0); %Door Mass
hObject.Data=RmParTableData;
function StopPush_CreateFcn(hObject, eventdata, handles)
hObject.UserData=0; %continue

%gui functions only for reset
function VarTable_CellEditCallback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)
function RmParTable_CellEditCallback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)
function NPntsEdit_Callback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)
function LelbowPop_Callback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)
function RelbowPop_Callback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)
function DirectionPop_Callback(hObject, eventdata, handles)
ResetRmPush_Callback(hObject, eventdata, handles)

%Dead gui functions
function Fig_CreateFcn(hObject, eventdata, handles)
function PauseTimeEdit_Callback(hObject, eventdata, handles)
function PauseTimeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PauseTimeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DirectionPop_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function NPntsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NPntsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function LelbowPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LelbowPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function RelbowPop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RelbowPop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MovieEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MovieEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MovieEdit_Callback(hObject, eventdata, handles)