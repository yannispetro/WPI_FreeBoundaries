n = 10;
A = -eye(n);
B = eye(n);
for i = 1:n-1
    A(i,i+1) = 1;
    B(i+1,i) = -1;
end

u = sym('u', [n,1]);
v = sym('v', [n,1]);

A*u
B*u
(A*u).^2.*(A*v) - (B*u).^2.*(B*v)
% D = (A*u).^2.*(A*v)