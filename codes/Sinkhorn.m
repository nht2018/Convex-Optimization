function [P, L_C] = Sinkhorn(C, a, b, epsilon, iterate_time)
    %   minimize    L_C(a,b) = <P,C> 
    %   subject to  P * 1_m == a
    %               P^T * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    %C:���۾���Ϊһ��n*m�ľ���
    %a:Ϊһ��n*1������,�Ҹ��������Ǹ�
    %b:Ϊһ��m*1������,�Ҹ��������Ǹ�
    %epsilon:�ͷ����ϵ��
    %iterate_time����������
    %����Ĳ���a,b������sum(a)==sum(b)
    [n, m] = size(C);
    
    %�ȱ�������ʹ��sum(a)==sum(b)==1
    lamda = sum(a);
    a = a / lamda;
    b = b / lamda;
    
    %��������Sinkhorn�㷨�����������:
    %   minimize    L_C^epsilon(a,b) = <P,C> - epsilon * H(P)
    %   subject to  P * 1_m == a
    %               P * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    % ����H(P) = -sum_{i,j} P_{ij}*log(P_{ij}-1) 
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
    
    %�ȱ������Ż�ԭ
    P = P * lamda;
    L_C = trace(P' * C);
end
    
    
    
    
    
    
    
    
    
    