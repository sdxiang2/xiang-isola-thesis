function [transition_degs] = Compare_Transition(t_fpd,c_fpd, p)
% Compares the transition degrees (degree at which half of the nodes are/
% are predicted to be isolas) between the theoretical distribution and
% continuted distribution.
% 
% Input:
%   t_fpd - (N-1)x4 matrix. Theoretical fraction-per-degree. The 
%       Summary_Stats variable created by erdos_renyi_distribution (we'll 
%       be using columns 1 and 3)
%   c_fpd - Maxdegree x 1 array, where Maxdegree is the maximum degree
%       observed in the randomly generated graphs. Continued 
%       fraction-per-degree. The isp_fpd variable output by 
%       ConStats_analysis2.
%
% Output:
%   transtion_degs - 1x2 array. The first element is the transition
%       degree of the theoretical distribution, and the second element is the
%       transition degree of the continued distribtuion.
%   Plots the two distributions together.


% ---------- BEGIN CODE ----------
% NOTE: Might need edits to work with BA graphs, especially the truncation
% part.

    
    n = size(t_fpd,1) + 1;

    % Finding transition degree for t_fpd
    below_point5 = t_fpd(:,3) < 0.5;
    t_transition_degree = find(below_point5,1,'last');
    leftpoint = [t_fpd(t_transition_degree,1), t_fpd(t_transition_degree,3) - 0.5]; % shifting down by 0.5 so we can find roots instead of intersection with y = 0.5
    rightpoint = [t_fpd(t_transition_degree+1,1), t_fpd(t_transition_degree+1,3) - 0.5];
    m =  (rightpoint(2) - leftpoint(2))/(rightpoint(1) - leftpoint(1));
    b = leftpoint(2) - m*leftpoint(1);
    transition_degs(1) = roots([m b])/n;

    % Finding transition degree for c_fpd
    below_point5 = c_fpd < 0.5;
    % Eliminating data points whose degree is in the top 10% -- the drop at
    % the end makes it harder to cleanly find the transition point
    cutoff_deg = icdf('Binomial',0.9,n,p);
    truncated_bp5 = below_point5(1:cutoff_deg);
    c_transition_degree = find(truncated_bp5,1,'last') - 1; % have to subtract 1 to get actual degree, since the array starts with degree 0
    leftpoint = [c_transition_degree, c_fpd(c_transition_degree+1) - 0.5];
    rightpoint = [c_transition_degree+1, c_fpd(c_transition_degree+2) - 0.5];
    m =  (rightpoint(2) - leftpoint(2))/(rightpoint(1) - leftpoint(1));
    b = leftpoint(2) - m*leftpoint(1);
    transition_degs(2) = roots([m b])/n;

    % Plotting t_fpd and c_fpd
    figure
    plot(t_fpd(:,1)/n, t_fpd(:,3), 'b-o');
    hold on
    plot((0:(length(c_fpd)-1))/n, c_fpd, 'r-*');
    legend('Theoretical', 'Continuation');
    xlabel('Degree/N');
    ylabel('Fraction isolas');





end