function [P, L_C] = log_domain_Sinkhorn(C, a, b, epsilon, iterate_time)
    %C:代价矩阵，为一个n*m的矩阵
    %a:为一个n*1的向量,且各分量均非负
    %b:为一个m*1的向量,且各分量均非负
    %epsilon:惩罚项的系数
    %iterate_time：迭代次数
    %输入的参数a,b需满足sum(a)==sum(b)
    [n, m] = size(C);
    
    %对偶问题为：
    %   maximize    L_C(a,b) = <f,a> + <g,b> 
    %   subject to  f_i + g_j <= C_{ij} for any i,j(1<=i<=n,1<=j<=m)            
    %f:为一个n*1的向量
    %g:为一个m*1的向量
       
    %等比例缩放使得sum(a)==sum(b)==1
    lamda = sum(a);
    a = a / lamda;
    b = b / lamda;
    
    K = exp( C /(-epsilon)); % K:Gibbs Kernal
    %接下来用log-domain-Sinkhorn算法解决如下问题,以避免Sinkhorn算法导致的
    %下溢。
    %   maximaize    L_C^epsilon(a,b) = <f,a> +<g,b> - epsilon * 
    %                                   <exp(f/epsilon),K * exp(g/epsilon)>
    %   subject to  P * 1_m == a
    %               P * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    
    function y = min_epsilon(z)
        z_ = min(z);
        y = z_-epsilon * log(sum(exp(-(z-z_)/epsilon)));
    end
    
    function y = Min_row_epsilon(A)
        y = zeros(size(A,1),1);
        for i = 1:size(A,1)
            y(i) = min_epsilon(A(i,:));
        end
    end
    
    function y = Min_col_epsilon(A)
        y = zeros(size(A,2),1);
        for i = 1:size(A,2)
            y(i) = min_epsilon(A(:,i));
        end
    end
    
    f = zeros(n,1);
    g = zeros(m,1);
    for i = 1 : iterate_time
        f = Min_row_epsilon(C -f -g') + f + epsilon * log(a);
        g = Min_col_epsilon(C - f - g') + g + epsilon * log(b);
    end
%     u_ = u .* min(a ./ (u .* (K * g)),1);
%     v_ = g .* min(b ./ (g .* (K' * u_)),1);
%     delta_a = a - u_ .* (K * v_);
%     delta_b = b - v_ .* (K' * u_);
%     P = diag(u_) * K * diag(v_) + delta_a * delta_b'/norm(delta_a,1);
    
    P = diag(exp(f/epsilon)) * K * diag(exp(g/epsilon));
    
    %等比例缩放还原
    P = P * lamda;
    L_C = trace(P' * C);
end
    
    
    
    
    
    
    
    
    
    