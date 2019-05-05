function [P, L_C] = log_domain_Sinkhorn(C, a, b, epsilon, iterate_time)
    %C:���۾���Ϊһ��n*m�ľ���
    %a:Ϊһ��n*1������,�Ҹ��������Ǹ�
    %b:Ϊһ��m*1������,�Ҹ��������Ǹ�
    %epsilon:�ͷ����ϵ��
    %iterate_time����������
    %����Ĳ���a,b������sum(a)==sum(b)
    [n, m] = size(C);
    
    %��ż����Ϊ��
    %   maximize    L_C(a,b) = <f,a> + <g,b> 
    %   subject to  f_i + g_j <= C_{ij} for any i,j(1<=i<=n,1<=j<=m)            
    %f:Ϊһ��n*1������
    %g:Ϊһ��m*1������
       
    %�ȱ�������ʹ��sum(a)==sum(b)==1
    lamda = sum(a);
    a = a / lamda;
    b = b / lamda;
    
    K = exp( C /(-epsilon)); % K:Gibbs Kernal
    %��������log-domain-Sinkhorn�㷨�����������,�Ա���Sinkhorn�㷨���µ�
    %���硣
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
    
    %�ȱ������Ż�ԭ
    P = P * lamda;
    L_C = trace(P' * C);
end
    
    
    
    
    
    
    
    
    
    