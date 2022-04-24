function [F,J,DFDP] = LDS_RHS(u,p,C)
% RHS with nonlinearity -μu + 2u^3 - u^5
%
% Inputs:
%   u - Nx1 array. Value of u at each node.
%   p - 2x1 or 1x2 array. Value of parameters. [d,μ]
%   C - NxN array. Graph laplacian.
%
% Outputs:
%   F - Nx1 array. RHS evaluated at u.
%   J - NxN array. Jacobian evaluated at u.
%   DFDP - Nx2 array. Derivative of function w.r.t. each parameter.


% ---------- BEGIN CODE ----------
  
    % rename parameters
    d = p(1);
    mu = p(2);
    N = size(u,1);

    % right-hand side
    F = d*C*u - mu*u + 2*u.^3 - u.^5;

    % Jacobian
    if nargout > 1
        J = d*C + spdiags(-mu + 6*u.^2 - 5*u.^4, 0, N, N);
    end

    % Derivative wrt parameters
    if nargout > 2
        DFDP = zeros(N,2);
        DFDP(:,1) = C*u; % derivative wrt d
        DFDP(:,2) = -u; % derivative wrt mu 
    end

end