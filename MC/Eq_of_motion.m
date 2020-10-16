function [ODEs] = Eq_of_motion(par,tf,ndof,targetPDF)

    ndof = par{1};

    m0   = par{2};
    c0   = par{3};
    k0   = par{4};

    e1   = par{5};
    e2   = par{6};

    S0   = par{7};

    D0    = par{8};

    M = m0*eye(ndof);
    C = c0*gallery('tridiag',ndof,-1,2,-1);
    K = k0*gallery('tridiag',ndof,-1,2,-1);

    Gk = e1*k0*eye(ndof);
    Gc = e2*c0*eye(ndof);

    D = 2*pi*S0*D0;
    
    x_str = cell(ndof,1);
    for i = 1:ndof
        x_str(i) = {['x' num2str(i) '(t)']};
    end

    syms t tf
    x(t) = str2sym(x_str);

    XD2 = M*diff(x,t,2) - C*diff(x,t,1) - K*x - Gc*diff(x,t,1).^3 - Gk*x.^3;

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
    Zd4 = matlabFunction(Zd4_symbolic,'Vars',{[X0;X1;X2;X3]});

    ODEs = @(z) [z(3:4);z(5:6);z(7:8);Zd4(z).'];

    free_BC = [2,3,4];

    % BCs
    BC_fixed = sym('BC_fixed', [2*ndof,1]);

    YB0 = sym('YB0', [ndof,1]);
    YB1 = sym('YB1', [ndof,1]);
    YB2 = sym('YB2', [ndof,1]);
    YB3 = sym('YB3', [ndof,1]);
    YB4 = sym('YB4', [ndof,1]);

    ind_tar = find(1:2*ndof == targetPDF);
    free_BC = [1:ind_tar-1 ind_tar+1:2*ndof];
    
    D_free = zeros(ndof,1);
    V_free = zeros(ndof,1);
    for j = 1:ndof
        if ismember(j,free_BC)
            D_free(j) = 1;
        end
        if ismember(j+ndof,free_BC)
            V_free(j) = 1;
        end
    end

    BCV = G2;
    BCD = G1 - diff(G2,t,1);

    BCV = subs( BCV, diff(x(t),t,4), YB4 );
    BCV = subs( BCV, diff(x(t),t,3), YB3 );
    BCV = subs( BCV, diff(x(t),t,2), YB2 );
    BCV = subs( BCV, diff(x(t),t,1), YB1 );
    BCV = subs( BCV, diff(x(t),t,0), YB0 );
    BCV_f = BCV(tf);

    BCD = subs( BCD, diff(x(t),t,4), YB4 );
    BCD = subs( BCD, diff(x(t),t,3), YB3 );
    BCD = subs( BCD, diff(x(t),t,2), YB2 );
    BCD = subs( BCD, diff(x(t),t,1), YB1 );
    BCD = subs( BCD, diff(x(t),t,0), YB0 );
    BCD_f = BCD(tf);
    
    % substitute the fixed (known) boundary values
    for j = 1:ndof
        if ~D_free(j)
            BCD_f = subs( BCD_f, YB0(j), BC_fixed(j) );
            BCV_f = subs( BCV_f, YB0(j), BC_fixed(j) );
        end
        if ~V_free(j)
            BCV_f = subs( BCV_f, YB1(j), BC_fixed(j+ndof) );
            BCD_f = subs( BCD_f, YB1(j), BC_fixed(j+ndof) );
        end
    end
    
    % For each ndof determine if free/fixed displacement and free/fixed
    % velocity and assign the proper boundary condition
    BBC = cell(2*ndof,1);
    for j = 1:ndof
        if D_free(j)
            BBC{j} = matlabFunction(BCD_f(j),'Vars',{[YB0;YB1;YB2;YB3;BC_fixed]});
        else
            BBC{j} = matlabFunction(YB0(j) - BC_fixed(j),'Vars',{[YB0;YB1;YB2;YB3;BC_fixed]});
        end
        
        if V_free(j)
            BBC{j+ndof} = matlabFunction(BCV_f(j),'Vars',{[YB0;YB1;YB2;YB3;BC_fixed]});
        else
            BBC{j+ndof} = matlabFunction(YB1(j) - BC_fixed(j+ndof),'Vars',{[YB0;YB1;YB2;YB3;BC_fixed]});
        end
    end
    BCs = @(y)cellfun(@(f)f(y),BBC);
    
end

% M = sym('M', [ndof, ndof]);
% C = sym('C', [ndof, ndof]);
% K = sym('K', [ndof, ndof]);
% 
% Gc = sym('Gc', [ndof, ndof]);
% Gk = sym('Gk', [ndof, ndof]);
% 
% M = diag(diag(M));
% Gc = diag(diag(Gc));
% Gk = diag(diag(Gk));