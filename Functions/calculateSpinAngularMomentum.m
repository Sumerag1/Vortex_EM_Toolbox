function sam = calculateSpinAngularMomentum(w0, mu, epsilon, Ex, Ey, Ez, Hx, Hy, Hz)
% 自动检测数据维度，并计算自旋角动量(SAM)
% 输入:
%   w0：频率
%   mu: 磁导率
%   epsilon: 介电常数
%   Ex, Ey, Ez: 电场分量
%   Hx, Hy, Hz: 磁场分量 (可选)
% 输出:
%   sam: 自旋角动量，保持与输入相同的维度

    % 检查输入参数数量
    if nargin < 5
        error('至少需要提供 mu, epsilon, Ex, Ey, Ez 五个参数');
    end
    % 检测数据维度
    originalSize = size(Ex);
    originalNdims = ndims(Ex);
    if originalNdims == 3
        fprintf('检测到3D数据，转换为列向量计算\n');
        % 将3D数据重塑为列向量
        Ex_vec = reshape(Ex, [], 1);
        Ey_vec = reshape(Ey, [], 1);
        Ez_vec = reshape(Ez, [], 1);
        
    % 如果有磁场分量，也进行重塑
        if nargin == 8
            Hx_vec = reshape(Hx, [], 1);
            Hy_vec = reshape(Hy, [], 1);
            Hz_vec = reshape(Hz, [], 1);
            % 计算自旋角动量
            sam_vec = calculateSAMCore(w0,mu, epsilon, Ex_vec, Ey_vec, Ez_vec, Hx_vec, Hy_vec, Hz_vec);
        else
            % 计算自旋角动量（只有电场）
            sam_vec = calculateSAMCore(w0,mu, epsilon, Ex_vec, Ey_vec, Ez_vec);
        end
        % 转换回原始形状
        sam = reshape(sam_vec, originalSize);
    elseif originalNdims == 1
        fprintf('检测到1D数据，直接计算\n');
        % 1D数据直接计算
        if nargin == 8
            sam = calculateSAMCore(mu, epsilon, Ex, Ey, Ez, Hx, Hy, Hz);
        else
            sam = calculateSAMCore(mu, epsilon, Ex, Ey, Ez);
        end
    end

end

% 核心函数
function sam_vec=calculateSAMCore(w0,epsilon,u,Ex,Ey,Ez,Hx,Hy,Hz)
% 自旋角动量核心计算函数
% 这里实现具体的SAM计算公式
    hasMagneticField = (nargin == 8);
    E = [Ex, Ey, Ez];
    E_conj = conj(E);
    % 检查是否有磁场分量
    if hasMagneticField==1
        H = [Hx, Hy, Hz];
        H_conj = conj(H);
        % 计算每个点的自旋密度S
        sam_vec = (1/(2*w0))*imag(epsilon.*cross(E_conj,E)+u.*cross(H_conj,H));
    else
        sam_vec = epsilon.*imag(cross(E_conj,E))./w0;
    end
end