function[fig] = Plot_BD(soln_struc,opt_struc)
% Plots and returns a bifurcation diagram figure.
%
% Inputs:
%   soln_struc - Solution structure output by TCLDS. Contains parameter and
%       solution data to be plotted.
%   opt_struc - (optional) Optional data structure output by TCLDS.
%       Contains fold, branch, and cutoff point data to be plotted.
%
% Outputs:
%   fig - figure containing bif. diagram.


% ---------- BEGIN CODE ----------

p = soln_struc.p;
sol = soln_struc.sol;
icp = soln_struc.icp;

sol = sol(any(sol,2),:);
p = p(any(sol,2),:);

plot(p(:,icp), vecnorm(sol'), 'LineWidth', 2);
ylabel('|u|');
xlabel('mu');
xlim([0,1]);

if nargin > 1
    fp = opt_struc.fold_pts;
    bp = opt_struc.branch_pts;
    cp = opt_struc.cutoff_pts;
    
    hold on
    
    for i = 1:size(fp,1)
        plot(p(fp(i),icp), norm(sol(fp(i),:)), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
    end
    
    for j = 1:size(bp,1)
        plot(p(bp(j),icp), norm(sol(bp(j),:)), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    end
    
    for k = 1:size(cp,1)
        plot(p(cp(k),icp), norm(sol(cp(k),:)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
    end
end

fig = gcf;

