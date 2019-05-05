function [P, L_C] = Sinkhorn(C, a, b, epsilon, iterate_time)
    %   minimize    L_C(a,b) = <P,C> 
    %   subject to  P * 1_m == a
    %               P^T * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    %C:代价矩阵，为一个n*m的矩阵
    %a:为一个n*1的向量,且各分量均非负
    %b:为一个m*1的向量,且各分量均非负
    %epsilon:惩罚项的系数
    %iterate_time：迭代次数
    %输入的参数a,b需满足sum(a)==sum(b)
    [n, m] = size(C);
    
    %等比例缩放使得sum(a)==sum(b)==1
    lamda = sum(a);
    a = a / lamda;
    b = b / lamda;
    
    %接下来用Sinkhorn算法解决如下问题:
    %   minimize    L_C^epsilon(a,b) = <P,C> - epsilon * H(P)
    %   subject to  P * 1_m == a
    %               P * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    % 其中H(P) = -sum_{i,j} P_{ij}*log(P_{ij}-1) 
    %          == -<P,log(P)-1>
    K = exp( C /(-epsilon)); % K:Gibbs Kernal
    v = ones(m,1);
    for i = 1 : iterate_time
        u = a ./ (K * v);
        v = b ./ (K' * u);
    end
    u_ = u .* min(a ./ (u .* (K * v)),1);
    v_ = v .* min(b ./ (v .* (K' * u_)),1);
    delta_a = a - u_ .* (K * v_);
    delta_b = b - v_ .* (K' * u_);
    P = diag(u_) * K * diag(v_) + delta_a * delta_b'/norm(delta_a,1);
    
    %等比例缩放还原
    P = P * lamda;
    L_C = trace(P' * C);
end
    
    
    
    
    
    
    
    
    
    