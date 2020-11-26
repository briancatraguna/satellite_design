function [sys,x0,str,ts] = saeroanim(t,x,u,flag,Config)
% animasi S-Function for displaying 6DoF attitude
% IRMUKa 8 Oktober 2006
% Copyright 1990-2003 The MathWorks, Inc.
% J.Hodgson  
% Modified by: S.Gage
% $Revision: 1.2.2.4 $  $Date: 2004/04/06 01:04:21 $
% Last modified w/o permission by: Harapan Rachman
% 8 Oktober 2006

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
     [sys,x0,str,ts]=mdlInitializeSizes(Config);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case {1 , 3, 9},
     sys=[];
     
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
     sys = [];
     if Config.Animenable
        mdlUpdate(t,x,u,Config);
     end

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u,Config);

otherwise
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%

   error('aeroblks:saeroanim:invalidflag',['Unhandled flag = ',num2str(flag)]);

end

% end sanim
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(Config)

%
% Set Sizes Matrix
%
sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 0;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialise the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialise the array of sample times
%
ts  = [-2 0]; % variable sample time

if ~Config.Animenable
   return
end

%
% Initialise Figure Window
%

   h_f=findobj('type','figure','Tag','6DOF anim');
   
   if isempty(h_f)
     h_anim=figure;
   else
     h_anim=h_f;
   end

   set(h_anim,'name','Animation Figure', ...
           'renderer','z','resize','off', ...
           'position',[253 112 546 449],'clipping','off', ...
           'Tag','6DOF anim');
   
   if ~isempty(h_anim)
      h_del = findobj(h_anim,'type','axes');
      delete(h_del);
      figure(h_anim);
   end
%
% Initialize Axes
%
   handle.axes(1)=axes;     % Satelit
   handle.axes(2)=axes;     % tak orbit
   handle.axes(3)=axes;     % tak benda
   set([handle.axes(1) handle.axes(2) handle.axes(3)],...
        'visible','off',...
        'dataaspectratio',[1 1 1],'projection','pers',...
        'units','normal', ...
        'position',[0.1 0.1 0.75 0.75], ...
        'Color',[0 0 0], ...
        'drawmode','fast', ...
        'clipping','off');
   b = 1.1*Config.craft;
   axis([-b b -b b -b b]);
   c = Config.craft;
   view([c c c]);

% ======================  Menggambar tak orbit  =====================
        Mo_transform = hgtransform;
       panah_satelit = hgtransform;
        num=5;		     % Number of z-y segments to make the circle
    %	count=1;	     % Counter for number of individual patches
        xtak = cell(19,1); %Preallocate for current number of patches
        ytak = cell(19,1); %Preallocate for current number of patches
        ztak = cell(19,1); %Preallocate for current number of patches

	theta=(360/2/num:360/num:(360+360/2/num))*pi/180;

	len 		= 2.5*Config.craft;       % Total Length (no units)
	radius      = len/100;              % Radius of body
 	s_fore      = len/10;				% Start of main body (w.r.t. nose)
   	thr_len 	= 0;                    % Length of Motor exit
   	l_fore      =len-s_fore-thr_len;	% Length of main body
   	c_g 		= len;                  % Position of c.g. w.r.t nose

    % ========================== sumbu x orbit ============================
    % Fore Body Shape
    %
	yc_range =  radius*sin(theta);
	zc_range = -radius*cos(theta);
	for i = 1:num
	   xtak{i}=[s_fore s_fore s_fore+l_fore s_fore+l_fore ]-c_g;
	   ytak{i}=[yc_range(i) yc_range(i+1) yc_range(i+1) yc_range(i)];
	   ztak{i}=[zc_range(i) zc_range(i+1) zc_range(i+1) zc_range(i)];
    end
    count=num+1;   
    %
    % Nose shape
    %
	for i = 1:num
	   xtak{count}=[s_fore s_fore 0 s_fore]-c_g;
	   ytak{count}=[yc_range(i) yc_range(i+1) 0 yc_range(i)];
      ztak{count}=[zc_range(i) zc_range(i+1) 0 zc_range(i)];
      count=count+1;
    end
    % Combine individual objects into a single set of co-ordinates
    % and roll through 45 degrees.
    %
    x=[];y=[];z=[];
    roll = [1         0                0;
            0   cos(45/180*pi)   sin(45/180*pi);
            0  -sin(45/180*pi)   cos(45/180*pi)];
    for i = 1:count
        x = [x xtak{i}'];
        y = [y ytak{i}'];
        z = [z ztak{i}'];
    end
   
    for i =1:4
        dum = [x(i,:);y(i,:);z(i,:)];
        dum = roll*dum;
        x(i,:)=-dum(1,:);
        y(i,:)=dum(2,:);
        z(i,:)=dum(3,:);
    end
    panahx = patch(x,y,z,[1 0 0]);
    panahx_satelit = patch(x,y,z,[1 0 0]);
% ========================== sumbu y orbit ============================
    % Fore Body Shape
    %
	xc_range =  radius*sin(theta);
	zc_range = -radius*cos(theta);
	for i = 1:num
	   xtak{i}=[xc_range(i) xc_range(i+1) xc_range(i+1) xc_range(i)];
	   ytak{i}=[s_fore s_fore s_fore+l_fore s_fore+l_fore ]-c_g;
	   ztak{i}=[zc_range(i) zc_range(i+1) zc_range(i+1) zc_range(i)];
    end
    count=num+1;
    %
    % Nose shape
    %
	for i = 1:num
	   xtak{count}=[xc_range(i) xc_range(i+1) 0 xc_range(i)];
	   ytak{count}=[s_fore s_fore 0 s_fore]-c_g;
       ztak{count}=[zc_range(i) zc_range(i+1) 0 zc_range(i)];
      count=count+1;
    end
    % =====================================================================
    % Combine individual objects into a single set of co-ordinates
    % and roll through 45 degrees.
    %
    x=[];y=[];z=[];
    pitch = [cos(45/180*pi)   0    -sin(45/180*pi);
             0                1     0;
             sin(45/180*pi)   0     cos(45/180*pi)];
    for i = 1:count
      x = [x xtak{i}'];
	  y = [y ytak{i}'];
	  z = [z ztak{i}'];
    end
   
    for i =1:4
      dum = [x(i,:);y(i,:);z(i,:)];
      dum = pitch*dum;
      x(i,:)=dum(1,:);
      y(i,:)=-dum(2,:);
      z(i,:)=dum(3,:);
    end
    panahy = patch(x,y,z,[0 1 0]);
    panahy_satelit = patch(x,y,z,[0 1 0]);

    % ========================== sumbu z orbit ============================
    % Fore Body Shape
    %
	xc_range =  radius*sin(theta);
	yc_range = -radius*cos(theta);
	for i = 1:num
	   xtak{i}=[xc_range(i) xc_range(i+1) xc_range(i+1) xc_range(i)];
	   ytak{i}=[yc_range(i) yc_range(i+1) yc_range(i+1) yc_range(i)];
	   ztak{i}=[s_fore s_fore s_fore+l_fore s_fore+l_fore ]-c_g;
    end
    count=num+1;   
    %
    % Nose shape
    %
	for i = 1:num
	   xtak{count}=[xc_range(i) xc_range(i+1) 0 xc_range(i)];
	   ytak{count}=[yc_range(i) yc_range(i+1) 0 yc_range(i)];
       ztak{count}=[s_fore s_fore 0 s_fore]-c_g;
      count=count+1;
    end
    % Combine individual objects into a single set of co-ordinates
    % and roll through 45 degrees.
    %
    x=[];y=[];z=[];
    yaw = [cos(45/180*pi)   sin(45/180*pi)   0;
          -sin(45/180*pi)   cos(45/180*pi)   0;
                0                 0          1];
    for i = 1:count
      x = [x xtak{i}'];
	  y = [y ytak{i}'];
	  z = [z ztak{i}'];
    end
   
    for i =1:4
      dum = [x(i,:);y(i,:);z(i,:)];
      dum = yaw*dum;
      x(i,:)=dum(1,:);
      y(i,:)=dum(2,:);
      z(i,:)=-dum(3,:);
    end
    panahz = patch(x,y,z,[0 0 1]);
    panahz_satelit = patch(x,y,z,[0 0 1]);
    
    set([panahx panahy panahz],'Parent',Mo_transform);
    set([panahx_satelit panahy_satelit panahz_satelit],...
         'Parent',panah_satelit);

    set([panahx_satelit panahy_satelit panahz_satelit
         panahx         panahy         panahz],'erasemode','nor',...
        'edgecolor',[0 0 0],'clipping','off');    
    
    set(handle.axes(2),'userdata',Mo_transform);
    set(handle.axes(3),'userdata',panah_satelit);
    
%   =============  selesai membuat tak orbit  =======================
%
% Draw in satellite shape
%
   [vertex_matrix, faces_matrix] = satellite_shape;
   handle.craft = patch('Vertices',vertex_matrix,...
          'Faces',faces_matrix,...
          'FaceVertexCData',hsv(6),...
          'FaceColor','flat');
   vert = Config.craft*get(handle.craft,'vertices');
   set(handle.axes(1),'userdata',vert)
   set(handle.craft,'erasemode','nor',...
       'edgecolor',[0 0 0],'clipping','off');
%
% Set Handles of graphics in Figure UserData
%
   set(h_anim,'userdata',handle);

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function mdlUpdate(t,x,u,Config)
%
% Obtain Handles of Figure Objects
%
    handle = get(findobj('type','figure','Tag','6DOF anim'),'userdata');

    if isempty(findobj('type','figure','Tag','6DOF anim'))
     %figure has been manually closed
     return
    end

%
% Form Transformation Matrix
%
    cph = cos(u(4));		% Roll
    sph = sin(u(4));
    cth = cos(u(5));		% Pitch
    sth = sin(u(5));
    cps = cos(u(6));		% Yaw
    sps = sin(u(6)); 
    attitude = [cth*cps sph*sth*cps-cph*sps cph*sth*cps+sph*sps
                cth*sps sph*sth*sps+cph*cps cph*sth*sps-sph*cps
                -sth         sph*cth         cph*cth];
%
% Update Craft Object 
%
   vert = get(handle.axes(1),'userdata');
   a=size(vert,1);
   dum =attitude*vert'+repmat(u(1:3,1),1,a);
   set(handle.craft,'vertices',dum');

 rx_angle = u(4);
 ry_angle = u(5);
 rz_angle = u(6);

 Rx = makehgtform('xrotate',rx_angle);
 Ry = makehgtform('yrotate',ry_angle);
 Rz = makehgtform('zrotate',rz_angle);

 Matrik_rotasi = Rz*Ry*Rx;    % Rotasi 3-2-1
 set(get(handle.axes(3),'userdata'),'Matrix',Matrik_rotasi);

%  ================  update tak orbit  =========================
%   Input sikap tak orbit
      G = 6.67E-11;
     Me = 5.98E24;
 tinggi = 400E3;
  Rbumi = 6371E3;
      a = tinggi + Rbumi;
      
       %wy = sqrt(G*Me/a^3);
       wy = 0;
  theta_y = wy*t;
    sudut = [0 theta_y 0]';

 rx_angle = sudut(1,1);
 ry_angle = sudut(2,1);
 rz_angle = sudut(3,1);

 Rx = makehgtform('xrotate',rx_angle);
 Ry = makehgtform('yrotate',ry_angle);
 Rz = makehgtform('zrotate',rz_angle);
 Sk = makehgtform('Scale',0.8);
 
 Matrik_rotasi = Rz*Ry*Rx*Sk;    % Rotasi 1-2-3,
 set(get(handle.axes(2),'userdata'),'Matrix',Matrik_rotasi);

%  =============================================================
%
% Set position of target view to Target
%
switch Config.camera_view

case 1,			% Fixed Observer Position
   set(handle.axes(1),'cameraupvector',[0 0 -1], ...
      'cameraposition',Config.camera_pos, ...
      'cameratarget',u(1:3), ...
      'cameraviewangle',Config.view);
case 2,			% Cockpit View
%   Target= Config.target';
   ax = Config.axes;
   seeker_dir = sqrt(sum((ax(2)-ax(1))^2+(ax(4)-ax(3))^2+(ax(6)-ax(5))^2))*attitude*[1;0;0];
   seeker_pos = attitude*[Config.craft 0 0]';
   set(handle.axes(1),'cameraupvector',attitude*[0 0 -1]', ...
           'cameraposition',u(1:3)+seeker_pos, ...
           'cameratarget',u(1:3)+seeker_dir, ...
           'cameraviewangle',Config.view);
case 3,			% Relative Position View
   set(handle.axes(1),'cameraupvector',[0 0 -1], ...
           'cameraposition',u(1:3)+Config.camera_pos', ...
           'cameratarget',u(1:3), ...
           'cameraviewangle',Config.view);
end

%
% Force MATLAB to Update Drawing
%
   drawnow

% end mdlUpdate

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u,Config)
    
    sys = ceil(t/Config.update)*Config.update+Config.update;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% satellite shape
% Function to draw shape of satellite object
%=============================================================================
%
function [vm,fm]=satellite_shape

% Membangun kubus
% Vertices       X   Y  Z
           vm = [-1 -1 -1;     % 1st vertex
                  1 -1 -1;     % ...
                  1  1 -1;
                 -1  1 -1;
                 -1 -1  1;
                  1 -1  1;
                  1  1  1;
                 -1  1  1];    % 8th vertex
% Faces
          fm = [1 2 6 5;    % 1st face
                2 3 7 6;    % ...
                3 4 8 7;
                4 1 5 8;
                1 2 3 4;
                5 6 7 8];   % 6th face
% ==============================================