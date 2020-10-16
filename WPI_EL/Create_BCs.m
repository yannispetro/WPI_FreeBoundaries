function BCs = Create_BCs(ndof,targetPDF,BC_free_sym)
    n_BCf = length(targetPDF);
    y = sym('y', [4*ndof,1]);
    BV = sym('BV', [n_BCf,1]);
    
    ind_fre = [];
    ind_fix = [];
    for j = 1:2*ndof
        if ismember(j,targetPDF)
            ind_fix = [ind_fix,j];
        else
            ind_fre = [ind_fre,j];
        end
    end
    
    for k = 1:n_BCf
        j = ind_fix(k);
        BC_free_sym = subs(BC_free_sym, y(j), BV(k));
    end
    
    BBC = cell(2*ndof,1);
    for j = ind_fre
        BBC{j} = matlabFunction(BC_free_sym(j),'Vars',{[y;BV]});
    end
    
    for k = 1:n_BCf
        j = ind_fix(k);
        BBC{j} = matlabFunction(y(j) - BV(k),'Vars',{[y;BV]});
    end
    
    BCs = @(y)cellfun(@(f)f(y),BBC);
    
    clear BBC
    
end

