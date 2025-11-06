function oam = calculateOrbitAngularMomentum(Ex, Ey, Ez, x, y, delta)
% 计算轨道角动量(OAM)密度
% 输入:
%   Ex, Ey, Ez: 电场分量（向量形式）
%   x, y: 坐标向量
%   delta: 网格间距
% 输出:
%   oam: 轨道角动量密度

% 重塑电场分量为网格格式
Nx = length(x);
Ny = length(y);
Ex_re=reshape(Ex,Nx,Ny)
%% derivative with x
z=fft(Ex_re);      % column fft
z=fftshift(z,1);   % column shift

n=0:Nx-1;
v=1i*pi/delta*(n-(Nx)/2)/((Nx)/2); % j*x in fourier space
z=z.*v.';  % column
z=ifftshift(z,1);  % column inverse shift

y_x=(ifft(z)); % column ifft


%% derivative with y
z=fft(Ex_re.');      % transpose and column fft
z=fftshift(z,1);     % column shift

n=0:Ny-1;
v=1i*pi/delta*(n-(Ny)/2)/((Ny)/2); % j*y in fourier space
z=z.*v.';  % column
z=ifftshift(z,1);  % column inverse shift
y_y=(ifft(z)); % column ifft

y_y=y_y.'; % transpose return back

%% ¶Ô PhiÇóµ¼

y_phi=y_x.*x-y_y.*y;
oam=conj(Ex_re).*(y_phi)*(-1i);
end
