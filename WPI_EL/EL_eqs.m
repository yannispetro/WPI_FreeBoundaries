function [ODEs,Lagrangian,BC_free_sym] = EL_eqs(par,tf)

    ndof = par.ndof;

    m0   = par.m0;
    c0   = par.c0;
    k0   = par.k0;
    
    w0   = par.w0;
    d0   = par.d0;

    e1   = par.e1;
    e2   = par.e2;

    Acos = par.Acos;
    wcos = par.wcos;

    S0   = par.S0;

    D0    = par.D0;

    M = m0*eye(ndof);
    C = c0*eye(ndof);
    K = w0*eye(ndof) + k0*full(gallery('tridiag',ndof,-1,2,-1));

    G1 = e1*eye(ndof);
    G2 = e2*eye(ndof);

    D = 2*pi*S0*(D0*D0.');
    
    x_str = cell(ndof,1);
    for i = 1:ndof
        x_str(i) = {['x' num2str(i) '(t)']};
    end

    syms t
    x(t) = str2sym(x_str);

    Fc0 = Acos*cos(wcos*t);
    Fc = [];
    for i = 1:ndof
        Fc = [Fc;Fc0];
    end
    
    % ODEs
    system = M*diff(x,t,2) + C*diff(x,t,1) + K*x + G1*(x./((x.^2 + d0^2).^(1.5))) - Fc;


    L = 1/2*system.'*inv(D)*system;

    X0 = sym('X0', [ndof,1]);
    X1 = sym('X1', [ndof,1]);
    X2 = sym('X2', [ndof,1]);
    X3 = sym('X3', [ndof,1]);
    X4 = sym('X4', [ndof,1]);

    LL = subs( L, diff(x(t),t,4), X4 );
    LL = subs( LL, diff(x(t),t,3), X3 );
    LL = subs( LL, diff(x(t),t,2), X2 );
    LL = subs( LL, diff(x(t),t,1), X1 );
    LL = subs( LL, diff(x(t),t,0), X0 );

    GG0 = gradient(LL, X0);
    GG1 = gradient(LL, X1);
    GG2 = gradient(LL, X2);

    G0 = subs(GG0, X4, diff(x(t),t,4));
    G0 = subs(G0, X3, diff(x(t),t,3));
    G0 = subs(G0, X2, diff(x(t),t,2));
    G0 = subs(G0, X1, diff(x(t),t,1));
    G0 = subs(G0, X0, diff(x(t),t,0));

    G1 = subs(GG1, X4, diff(x(t),t,4));
    G1 = subs(G1, X3, diff(x(t),t,3));
    G1 = subs(G1, X2, diff(x(t),t,2));
    G1 = subs(G1, X1, diff(x(t),t,1));
    G1 = subs(G1, X0, diff(x(t),t,0));

    G2 = subs(GG2, X4, diff(x(t),t,4));
    G2 = subs(G2, X3, diff(x(t),t,3));
    G2 = subs(G2, X2, diff(x(t),t,2));
    G2 = subs(G2, X1, diff(x(t),t,1));
    G2 = subs(G2, X0, diff(x(t),t,0));

    EL = G0 - diff(G1,t) + diff(G2,t,2);

    ELL = subs( EL, diff(x(t),t,4), X4 );
    ELL = subs( ELL, diff(x(t),t,3), X3 );
    ELL = subs( ELL, diff(x(t),t,2), X2 );
    ELL = subs( ELL, diff(x(t),t,1), X1 );
    ELL = subs( ELL, diff(x(t),t,0), X0 );

    Sol = solve( ELL == 0, X4 );
    V = fieldnames(Sol);

    Zd4_symbolic = [];
    for i = 1:ndof
        Zd4_symbolic = [Zd4_symbolic, getfield(Sol,V{i})];
    end
    Zd4 = matlabFunction(Zd4_symbolic,'Vars',{[X0;X1;X2;X3],t});

    ODEs = @(z,t) [z(ndof+1:4*ndof);Zd4(z,t).'];

    % Lagrangian
    Lagrangian = matlabFunction(LL,'Vars',{[X0;X1;X2;X3;t]});
    
    % BCs
    YB0 = sym('YB0', [ndof,1]);
    YB1 = sym('YB1', [ndof,1]);
    YB2 = sym('YB2', [ndof,1]);
    YB3 = sym('YB3', [ndof,1]);
    YB4 = sym('YB4', [ndof,1]);

    BCV = G2;
    BCD = G1 - diff(G2,t,1);

%     BCV = subs( BCV, diff(x(t),t,4), YB4 );
    BCV = subs( BCV, diff(x(t),t,3), YB3 );
    BCV = subs( BCV, diff(x(t),t,2), YB2 );
    BCV = subs( BCV, diff(x(t),t,1), YB1 );
    BCV = subs( BCV, diff(x(t),t,0), YB0 );
    BCV_f = BCV(tf);

%     BCD = subs( BCD, diff(x(t),t,4), YB4 );
    BCD = subs( BCD, diff(x(t),t,3), YB3 );
    BCD = subs( BCD, diff(x(t),t,2), YB2 );
    BCD = subs( BCD, diff(x(t),t,1), YB1 );
    BCD = subs( BCD, diff(x(t),t,0), YB0 );
    BCD_f = BCD(tf);
    
    y = sym('y', [4*ndof,1]);
    BC_free_sym = [subs(BCD_f,[YB0;YB1;YB2;YB3],y);subs(BCV_f,[YB0;YB1;YB2;YB3],y)];
    
end