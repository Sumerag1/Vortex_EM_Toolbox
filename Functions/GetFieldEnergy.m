function fieldEnergy=GetFieldEnergy(mu,Ex,Ey,Ez,Hx,Hy,Hz)
%得到场能量
% 输入:
%   w0：频率
%   mu: 磁导率
%   epsilon: 介电常数
%   Ex, Ey, Ez: 电场分量
%   Hx, Hy, Hz: 磁场分量 (可选)
% 输出:
%   fieldEnergy: 场能量
E=real(Ex.*conj(Ex)+Ey.*conj(Ey)+Ez.*conj(Ez));
B=real(mu.*(Hx.*conj(Hx)+Hy.*conj(Hy)+Hz.*conj(Hz)));
fieldEnergy=(E+B)/(8*pi);
end