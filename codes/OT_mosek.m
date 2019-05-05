function  [P, L_C] = OT_mosek(C, a, b)
 %   minimize    L_C(a,b) = <P,C> 
 %   subject to  P * 1_m == a
 %               P^T * 1_n == b
 %               P_{ij} >= 0 for any i,j (1<=i<=n,1<=j<=m)
[n, m] = size(C);
C = reshape(C',[n*m, 1]);
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
A = [U;V];
blc = [a;b];
buc = [a;b];
blx = zeros(m*n, 1);
bux = [];
[res] = msklpopt(C,A,blc,buc,blx,bux,[],'minimize');
P = res.sol.itr.xx;
L_C = P' * C;
P = reshape(P,[m, n])';
