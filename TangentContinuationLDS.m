function [soln_struc,opt_struc] = TangentContinuationLDS(L,fhandle,u0,pars,icp,ds,nmx,options)
% Uses tangent numerical continuation to calculate bifurcation diagrams.
%
% Inputs:
%   L - matrix. Graph laplacian.
%   fhandle - function handle. RHS of diff eq. Takes u (array, value of u 
%       at each node in graph) and pars (array, value of parameters) and 
%       returns F (array, value of function), J (matrix, Jacobian), and 
%       DFDP (matrix, derivative of function w.r.t. each parameter)
%   u0 - array. initial values of u at each node.
%   pars - array. initial parameter values.
%   icp - int. index of parameter in pars to be continued.
%   ds - number. step size. set to negative if continuing towards the left
%       first.
%   nmx - int. number of steps.
%   options - array. Contains the following options for continuation:
%       detect_branch - {0,1,2}. Switches branch point detection option. 
%           0: no detection. 1: fold point detection. 2: branch point
%           detection.
%       cutoff_val - number. “cutoff value” for parameter (the indices at 
%           which the parameter crosses this value are recorded, but the 
%           continuation will run to nmx). If cutoff_val is set, the
%           continuation will stop once 5 cutoff points are found. This
%           is to aid in figure 8 isola detection.
%
% Outputs:
%   soln_struc - array. [parameters, solutions, icp, L]. Contains main
%       continuation data.
%   opt_struc - array. [fold points, branch points, cutoff points].
%       Contains optional continuation data.


% ---------- BEGIN CODE ----------

  % Default detect_branch: 0. Default cutoff_val: -1.
  if nargin > 7
      db = options(1);
      cv = options(2);
  else
      db = 0;    
      cv = -1;
  end

  
  % Initialise variables
  ndim = size(u0,1);
  p = zeros(nmx+1,2); % parameter values at each step
  sol = zeros(nmx+1,ndim); % solution at each step
  v1 = zeros(ndim+1,1);
  ds1 = ds;
  fold_pts = zeros(0,1);
  branch_pts = zeros(0,1);
  cutoff_pts = zeros(0,1);
  
  
  % Options to the nonlinear solver
  opts = optimset('Display','off','Jacobian','on');
  opts2 = optimset('Display','off','Algorithm','levenberg-marquardt'); % this is for fold detection's fsolve

  
  % Converge initial guess
  v1(ndim+1) = pars(icp);
  v1(1:ndim) = fsolve( @(u) fhandle(u,pars), u0, opts);
  p(1,:) = pars;  
  sol(1,:) = v1(1:ndim)';

  
  % Predictor
  pars(icp) = v1(ndim+1);
  [F,DFDU,DFDP] = fhandle(v1(1:ndim),pars);
  DFDP = DFDP(:,icp); 
  sec = null(full([DFDU, DFDP]));
  sec = sign(sec(ndim+1))*sign(ds)*sec;
  v = v1 + abs(ds1)*sec;
  
  
  % Calculate initial tangent vector (for fold pt detection) and initial
  % determinant of Jacobian(for branch pt detection)  
  F_U = [DFDU,DFDP];
  W = fsolve(@(W) F_U * W, sec, opts2);
  W_series = zeros(3,1);
  W_series(1) = W(ndim+1);
  D_old = det(DFDU);  

  
  % Start tangent continuation
  for n = 1:nmx
      
    % Corrector
    [v,~,exitflag] = fsolve( @(v) SecantCorrector(v), v, opts);
	if exitflag<=0
		fprintf('%2d %14.4e\n', exitflag, bd(n+1,2));
    end
    
    
    % Bookkeeping
    pars(icp) = v(ndim+1);
    p(n+1,:) = pars;
    sol(n+1,:) = v(1:ndim)';
	v1 = v;
        
        
    % Detecting when parameter crosses cutoff
    if (cv - p(n,icp)) * (cv - p(n+1,icp)) <= 0
        % finding exact cutoff point
        v(ndim+1) = cv;
        pars(icp) = cv;
        v(1:ndim) = fsolve(@(u) fhandle(u(1:ndim), pars), v(1:ndim));
        % Bookkeeping
        p(n,:) = pars;
        sol(n,:) = v(1:ndim)';
        cutoff_pts = [cutoff_pts; n];
        if size(cutoff_pts,1) == 5
            break
        end
    end    

    
	% Predictor
	[F,DFDU,DFDP] = fhandle(v1(1:ndim),pars);
    DFDP = DFDP(:,icp);
    sec1 = null(full([DFDU, DFDP]));
	sec = sign(sec'*sec1)*sec1;
	v = v1 + abs(ds1)*sec;
    
    
    % Detecting fold point/branch point
    switch db
        case 1
            F_U = [DFDU,DFDP];
            W = fsolve(@(W) F_U * W, sec, opts2);
            switch n
                case 1
                    W_series(2) = W(ndim+1);
                case 2
                    W_series(3) = W(ndim+1);
                otherwise
                    %if W_series(1) * W_series(2)^2 * W_series(3) <= 0
                    if W_series(2) * W_series(3) <= 0  
                        % Finding exact fold point
                        v1 = fsolve(@(v) bp_func(v), v1);
                        % Bookkeeping
                        p(n,icp) = v1(ndim+1);
                        sol(n,:) = v1(1:ndim);
                        fold_pts = [fold_pts;n];
                    end
                    W_series(1:2) = W_series(2:3);
                    W_series(3) = W(ndim+1);
            end
        case 2
            D_new = det(DFDU);
            if D_new * D_old <= 0
                % Finding exact fold point
                v1 = fsolve(@(v) bp_func(v), v1);
                % Bookkeeping
                p(n,icp) = v1(ndim+1);
                sol(n,:) = v1(1:ndim);
                branch_pts = [branch_pts;n];
            end
            D_old = D_new;
    end       

  end
  
  
  % Packaging into return structures
  soln_struc.p = p;
  soln_struc.sol = sol;
  soln_struc.icp = icp;
  soln_struc.L = L;
  opt_struc.fold_pts = fold_pts;
  opt_struc.branch_pts = branch_pts;
  opt_struc.cutoff_pts = cutoff_pts;
  
  
  function [G,DG] = SecantCorrector(v)
    % Compute F
    pars(icp) = v(ndim+1);
	F = fhandle(v(1:ndim),pars);

    % Extended system
    G = [F; sec'*(v-v1) - abs(ds1)];
    
    if nargout > 1
      % Jacobian of F
      [F,DFDU,DFDP] = fhandle(v(1:ndim),pars);
      
      DFDP = DFDP(:,icp);

      % Jacobian of the extended system
      DG = [DFDU, DFDP; sec'];
    end
  end

  
    function[H] = bp_func(v)
        pars(icp) = v(ndim+1);
        [F,J] = fhandle(v(1:ndim),pars);
        H = [F; det(J)];
    end

end
