function [patterns] = GetPatterns(soln_struc, opt_struc)
% Returns pattern data at cutoff points.
%
% Inputs:
%   soln_struc - Solution structure output by TCLDS.
%   opt_struc - Optional data structure output by TCLDS.
%
% Outputs:
%   patterns - pattern data at cutoff points. Formatted as
%       [u1, ..., uN, cont. parameter]


% ---------- BEGIN CODE ----------

    cp = opt_struc.cutoff_pts;
    p = soln_struc.p;
    sol = soln_struc.sol;
    icp = soln_struc.icp;
    
    ndim = size(sol,2);
    npat = size(cp,1);
    patterns = zeros(npat,ndim+1);
    
    for i = 1:npat
        patterns(i,1:ndim) = sol(cp(i),:);
        patterns(i,ndim+1) = p(cp(i),icp);
    end

end