function [P, L_C] = OT_gurobi(C, a, b)
    %   minimize    L_C(a,b) = <P,C> 
    %   subject to  P * 1_m == a
    %               P^T * 1_n == b
    %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
    %C:���۾���Ϊһ��n*m�ľ���
    %a:Ϊһ��n*1������,�Ҹ��������Ǹ�
    %b:Ϊһ��m*1������,�Ҹ��������Ǹ�
    %����Ĳ���a,b������sum(a)==sum(b)
    [n, m] = size(C); 
    model.modelsense = 'Min';
    model.modelname  ='OT_gurobi';
    c = reshape(C',[1, n*m]);
    model.obj = c;
    
    U = zeros(n, n*m);
    for i = 1 : n
        for j = 1 : m
            U(i,m*(i-1)+j) = 1; 
        end
    end
    V = [];
    for i = 1:n
        V = [V, eye(m)];
    end
    model.A = sparse([U; V]);
    model.rhs = [a; b];
    model.lb = zeros(n*m,1);
    sense = [];
    for i = 1:(n+m)
        sense = [sense, '='];
    end
    model.sense = sense;
    result = gurobi(model);
     
   	P = reshape(result.x,[m, n])';
    L_C = trace(P' * C);   
end
