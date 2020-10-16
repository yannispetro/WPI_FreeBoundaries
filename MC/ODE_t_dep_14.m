function z_dot = ODE_t_dep_14(t,z,par,Ce,Ke,tt,At1,At2)

    m1 = par(1);
    m2 = par(2);
    c1 = par(3);
    c2 = par(4);
    k1 = par(5);
    k2 = par(6);

    ec1 = par(7);
    ec2 = par(8);
    ek1 = par(9);
    ek2 = par(10);
    
    S0 = par(11);

    M = [m1 0; 0 m2];
    C = [c1+c2 -c2; -c2 c2];
    K = [k1+k2 -k2; -k2 k2];

    Cet = zeros([2,2]);
    Cet(1,1) = interp1(tt,squeeze(Ce(1,1,:)),t,'linear');
    Cet(1,2) = interp1(tt,squeeze(Ce(1,2,:)),t,'linear');
    Cet(2,1) = interp1(tt,squeeze(Ce(2,1,:)),t,'linear');
    Cet(2,2) = interp1(tt,squeeze(Ce(2,2,:)),t,'linear');

    Ket = zeros([2,2]);
    Ket(1,1) = interp1(tt,squeeze(Ke(1,1,:)),t,'linear');
    Ket(1,2) = interp1(tt,squeeze(Ke(1,2,:)),t,'linear');
    Ket(2,1) = interp1(tt,squeeze(Ke(2,1,:)),t,'linear');
    Ket(2,2) = interp1(tt,squeeze(Ke(2,2,:)),t,'linear');
    
    A1 = interp1(tt,At1,t,'linear');
    A2 = interp1(tt,At2,t,'linear');
    
    A = [A1,0;0,A2];
    
    nt = length(tt);

    Minv = inv(M);
    
    z_dot = zeros([14,1]);

    z_dot(1) = z(3);
    z_dot(2) = z(4);
    
    g1 = c1*z(3)^3*ec1+ek1*k1*z(1)^3+5*ek1*k1*z(1)*z(5)+5*c1*z(3)*ec1*z(12);
    g2 = c2*z(4)^3*ec2+ek2*k2*z(2)^3+5*ek2*k2*z(2)*z(9)+5*c2*z(4)*ec2*z(14);
    
    temp1 = - Minv*C*[z(3);z(4)] - Minv*K*[z(1);z(2)] - Minv*[g1;g2];
    z_dot(3) = temp1(1);
    z_dot(4) = temp1(2);
    
    V = Vec10toMat4(z(5:end));
    
    Gt = [zeros([2,2]), eye(2); -Minv*(K+Ket), -Minv*(C+Cet)];
    D = 2*pi*S0*eye(2);
    Lt = A*D*A.';
    Th = [zeros([2,2]), zeros([2,2]); zeros([2,2]), Minv*Lt*transpose(Minv)];
    
    V_dot = Gt*transpose(V) + V*transpose(Gt) + Th;
    
    z_dot(5:end) = Mat4toVec10(V_dot);
    
end

%Reshaping test
% aa = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
% bb = reshape(aa,[4 4]);
% cc = reshape(bb,[16 1]);
