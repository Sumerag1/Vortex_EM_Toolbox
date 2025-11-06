function sam = calculateSpinAngularMomentum(w0, epsilon, mu, Ex, Ey, Ez, Hx, Hy, Hz)
% 自动检测数据维度，并计算自旋角动量(SAM)
% 输入:
%   w0：频率
%   epsilon: 介电常数  
%   mu: 磁导率
%   Ex, Ey, Ez: 电场分量
%   Hx, Hy, Hz: 磁场分量 (可选)
% 输出:
%   sam: 自旋角动量，保持与输入相同的维度

    % 检查输入参数数量
    if nargin < 6
        error('至少需要提供 w0, epsilon, mu, Ex, Ey, Ez 6个参数');
    end
    
    % 检测数据维度
    originalSize = size(Ex);
    originalNdims = ndims(Ex);
    
    fprintf('数据维度: %s, 维数: %d\n', mat2str(originalSize), originalNdims);
    
    % 初始化sam（重要！确保所有路径都有返回值）
    sam = [];
    
    if originalNdims == 3
        fprintf('检测到3D数据，转换为列向量计算\n');
        % 将3D数据重塑为列向量
        Ex_vec = reshape(Ex, [], 1);
        Ey_vec = reshape(Ey, [], 1);
        Ez_vec = reshape(Ez, [], 1);
        
        % 如果有磁场分量，也进行重塑
        if nargin == 9
            Hx_vec = reshape(Hx, [], 1);
            Hy_vec = reshape(Hy, [], 1);
            Hz_vec = reshape(Hz, [], 1);
            sam_vec = calculateSAMCore(w0, epsilon, mu, Ex_vec, Ey_vec, Ez_vec, Hx_vec, Hy_vec, Hz_vec);
        else
            sam_vec = calculateSAMCore(w0, epsilon, mu, Ex_vec, Ey_vec, Ez_vec);
        end
        % 转换回原始形状
        sam = reshape(sam_vec, originalSize);
        
    elseif originalNdims == 2
        % 处理2D数据（包括列向量和矩阵）
        fprintf('检测到2D数据\n');
        if iscolumn(Ex)  % 如果是列向量 [N, 1]
            fprintf('数据是列向量，直接计算\n');
            if nargin == 9
                sam = calculateSAMCore(w0, epsilon, mu, Ex, Ey, Ez, Hx, Hy, Hz);
            else
                sam = calculateSAMCore(w0, epsilon, mu, Ex, Ey, Ez);
            end
        else  % 如果是2D矩阵 [M, N]
            fprintf('数据是2D矩阵，转换为列向量计算\n');
            Ex_vec = reshape(Ex, [], 1);
            Ey_vec = reshape(Ey, [], 1);
            Ez_vec = reshape(Ez, [], 1);
            
            if nargin == 9
                Hx_vec = reshape(Hx, [], 1);
                Hy_vec = reshape(Hy, [], 1);
                Hz_vec = reshape(Hz, [], 1);
                sam_vec = calculateSAMCore(w0, epsilon, mu, Ex_vec, Ey_vec, Ez_vec, Hx_vec, Hy_vec, Hz_vec);
            else
                sam_vec = calculateSAMCore(w0, epsilon, mu, Ex_vec, Ey_vec, Ez_vec);
            end
            sam = reshape(sam_vec, originalSize);
        end
        
    elseif originalNdims == 1
        fprintf('检测到1D数据，直接计算\n');
        if nargin == 9
            sam = calculateSAMCore(w0, epsilon, mu, Ex, Ey, Ez, Hx, Hy, Hz);
        else
            sam = calculateSAMCore(w0, epsilon, mu, Ex, Ey, Ez);
        end
    else
        error('不支持的维度: %d', originalNdims);
    end
    
    % 最终检查
    if isempty(sam)
        error('计算失败：sam为空，请检查数据维度处理逻辑');
    end
    
    fprintf('计算结果维度: %s\n', mat2str(size(sam)));
end

% 核心函数保持不变
function sam_vec = calculateSAMCore(w0, epsilon, mu, Ex, Ey, Ez, Hx, Hy, Hz)
% 自旋角动量核心计算函数
    hasMagneticField = (nargin == 9);
    
    % 确保是列向量
    Ex = Ex(:); Ey = Ey(:); Ez = Ez(:);
    
    E = [Ex, Ey, Ez];
    E_conj = conj(E);
    
    % 检查是否有磁场分量
    if hasMagneticField
        Hx = Hx(:); Hy = Hy(:); Hz = Hz(:);
        H = [Hx, Hy, Hz];
        H_conj = conj(H);
        % 计算每个点的自旋密度S
        sam_vec = (1/(2*w0)) * imag(epsilon .* cross(E_conj, E) + mu .* cross(H_conj, H));
    else
        % 只有电场的自旋角动量
        sam_vec = (epsilon/(2*w0)) .* imag(cross(E_conj, E));
    end
end
