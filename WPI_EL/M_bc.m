function res = M_bc(ya, yb, ODEs, BCs, BVX)

    nsdof = length(yb)/2;

    res = zeros([1,2*nsdof]);
    for i = 1:nsdof
        res(i) = ya(i) - 0;
    end

    res(nsdof+1:end) = BCs([yb;BVX]).';

end