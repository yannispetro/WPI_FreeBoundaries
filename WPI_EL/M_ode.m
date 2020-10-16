function z_dot= M_ode(t, z, ODEs, BCs, BVX)

    z_dot = ODEs(z,t);

end