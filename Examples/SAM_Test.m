
clc;clear

%% 参数设置
wavelength=1;  %  ***wavelength
k0=2*pi/wavelength;  %  wavenumber
c=3e8;
w0=c/wavelength*2*pi;
Z=120*pi;  % wave impedance
Volume=1;
%an=(Volume/(pi*4/3)).^(1/3);
TL=2;   % Topological charge   ***

%% number and position of source points
N=12;   %  *** number of point sources

% radius and angle
radius= 0.5;  %  *** radius
theta=linspace(0,2*pi-2*pi/(N),N);
[x00,y00]=pol2cart(theta,radius);

%% sampling the observation area xoy plane
%  z direction slices
delta=0.4;       % step ***
size=30;         % xoy range ***

% xoy plane
xo=-size:delta:size-delta;
yo=-size:delta:size-delta;
zo=1;  % observation z  ***

% size
Nx=length(xo);
Ny=length(yo);

%  position of observation point
index=0;
for m=1:Nx;  %  X
    for n=1:Ny;  %  Y
            index=index+1;
            posx(index)=xo(m);
            posy(index)=yo(n);
            posz(index)=zo(1);
    end
end
N_p=index;

figure(1);
plot(x00,y00,'r.')
hold on
plot([x00,0.5],[y00,0],'b--')
xlabel('x (a.u.)')
ylabel('y (a.u.)')
axis equal
% N_p=index;
title('Location of point sources')

%% dyadic Green function
Greenxx=zeros(N_p,N);
Greenyy=zeros(N_p,N);
Greenzz=zeros(N_p,N);
Greenxy=zeros(N_p,N);
Greenxz=zeros(N_p,N);
Greenyz=zeros(N_p,N);

for index_f=1:N_p; % 场点个数
    %  计数
    %  N_p-index_f
    
    %  观测点位置
    x=posx(index_f);
    y=posy(index_f);
    z=posz(index_f);
    
    for index_s=1:N;
        
        %  源点位置
        xx=x00(index_s);
        yy=y00(index_s);
        zz=0;
        
        %  距离
        R=sqrt((xx-x)^2+(yy-y)^2+(zz-z)^2);
        alpha=k0*R;
        %  方向
        cosx=(x-xx)/R;
        cosy=(y-yy)/R;
        cosz=(z-zz)/R;
        %  常数
        const1=-j*k0*Z*k0*Volume*exp(-j*alpha)/(4*pi*alpha^3);
        const2=3-alpha^2+3*j*alpha;
        const3=(alpha)^2-1-j*alpha;
        Greenxx(index_f,index_s)=const1*(const3+cosx*cosx*const2);
        Greenyy(index_f,index_s)=const1*(const3+cosy*cosy*const2);
        Greenzz(index_f,index_s)=const1*(const3+cosz*cosz*const2);
        Greenxy(index_f,index_s)=const1*cosx*cosy*const2;
        Greenxz(index_f,index_s)=const1*cosx*cosz*const2;
        Greenyz(index_f,index_s)=const1*cosy*cosz*const2;
        
    end
end

%  并矢格林函数
G=[Greenxx Greenxy Greenxz;
    Greenxy Greenyy Greenyz;
    Greenxz Greenyz Greenzz];

%  电流源
for index_s=1:N;
    Jx(index_s)=exp(i*TL*theta(index_s));  % 只在x方向有电流源
    Jy(index_s)=0;
    Jz(index_s)=0;
end

J=[Jx Jy Jz].';

% 场的结果
Etot=G*J;

Ex=Etot([1:N_p],1);
Ey=Etot([N_p+1:2*N_p],1);
Ez=Etot([2*N_p+1:3*N_p],1);

SAM=calculateSpinAngularMomentum(w0,epsilon,u,Ex,Ey,Ez);

imagesc(SAM(:,:,3));%得到的是二维平面的数据，此处展示Sz的结果
