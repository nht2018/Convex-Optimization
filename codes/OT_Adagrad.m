function varargout = OT_Adagrad(C, mu, nu, rou, step_size, iterate_time,tol)
    TOL = tol;
    [m, n] = size(C);
    %   minimize    <Pi,C> 
    %   subject to  Pi * 1_ == mu
    %               Pi^T * 1_n == nu
    %               Pi_{ij} >= 0 for any i,j (1<=i<=m,1<=j<=n)
    %C:���۾���Ϊһ��n*m�ľ���
    %mu:Ϊһ��m*1������,�Ҹ��������Ǹ�
    %nu:Ϊһ��n*1������,�Ҹ��������Ǹ�
    %rou:�ͷ����ϵ��
    %iterate_time����������
    %����Ĳ���mu,nu������sum(mu)==sum(nu)
   
    
    %�����������պ���
    %L(Pi,y,z)=<Pi,C>-y'*(Pi*ones(n,1)-mu)-z'*(Pi'*ones(m,1)-nu)
    %+0.5*rou*norm(Pi*ones(n,1)-mu,2)+0.5*rou*norm(Pi'*ones(m,1)-nu,2)
    
    %ADMM����
    gamma = 1;
    y = zeros(m,1);
    z = zeros(n,1);
    Pi = rand(m,n);
    i = 1;
    value = trace(Pi' * C);
    r = 0;
    mini = 1e-7;
    dec = -1;
    while i < iterate_time && ~(abs(dec) < value * TOL)
        nabla_L = C - y - z' + rou* (sum(Pi,1) + sum(Pi,2)- mu - nu');
        r = r+ sum(sum(nabla_L.^2));
        Pi = max(Pi - step_size/(mini + r^0.5)*nabla_L ,0);%��һ���ݶȲ�
        y = y - gamma * rou * (Pi * ones(n,1) - mu);
        z = z - gamma * rou * (Pi' * ones(m,1) - nu);
        dec = value - trace(Pi' * C);
        value = value - dec;
        if nargout == 3
            out(i) = value;
        end
        i = i+ 1;
    end
    %�ȱ������Ż�ԭ
    MinCost = trace(Pi' * C);
    varargout{1} = Pi;
    varargout{2} = MinCost;
    if nargout == 3
        varargout{3} = out;
    end
end
    
    
    
    
    
    
    
    
    
    