function z_dot = ODE_MC( t,z, My,C,K,d0,G1,G2,tt,F )
    
    s = size(F);
    N = s(1);   
    
    Z1 = z(1:N);
    Z2 = z(N+1:end);
 
    Ft = interp1(tt,F.',t,'linear').';

    z_dot = zeros(size(z));
    z_dot(1:N) = Z2;
    z_dot(N+1:end) = My*Ft - My*( C*Z2 + K*Z1 + G1*(Z1./((Z1.^2 + d0^2).^(1.5))) );
%     z_dot(N+1:end) = My*Ft - My*(C*Z2 + K*Z1 + G1*(Z1.^3) - G2*( (GA*Z1).^2.*(GA*Z2) - (GB*Z1).^2.*(GB*Z2) ) );

end

