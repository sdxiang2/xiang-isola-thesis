function [] = DispPatterns(soln_struc, patterns)
% Interactively displays pattern data.
%
% Inputs
%   soln_struc - Solution structure output by TCLDS. Used to plot bifurca-
%       tion diagram.
%   patterns - Output of GetPatterns. Contains pattern data corrresponding
%       to soln_struc.
% Outputs
%   None. Creates figure window with graph, bif. diagram, and grid. Input
%   1 to the command line to go forwards, 0 to go backwards, and 'q' to
%   quit.


% ---------- BEGIN CODE ----------
    
    % Creating graph object from L
    L = soln_struc.L;
    adj_matrix = -(L + diag(-diag(L)));
    G = graph(adj_matrix);


    % Number of nodes
    N = size(L,1);
    
    
    % Number of intersections
    M = size(patterns,1);
    
    
    % Pattern display loop
    j = 0;
    key = 1; % 1 to go forwards, 0 to go backwards
    figure('Position', [300 300 700 700]);
    while key ~= 'q'
        if key == 1 && j < M
            j = j + 1;
        elseif key == 1 && j == M
            j = 1;
        elseif key == 0 && j > 1
            j = j - 1;
        elseif key == 0 && j == 1
            j = M;
        end
        
        % data is what's displayed at each new key input
        data = ones(M,N);
        for r = 1:j
            for c = 1:N % loops over u1 ... uN
                if patterns(r,c) > 0.8 % if u+
                    data(r,c) = 3;
                elseif patterns(r,c) > 0.3 % if u-
                    data(r,c) = 2;
                end
            end
        end
        
        % Identifying u+, u- nodes to highlight on the graph
        red = [];
        yellow = [];
        for i = 1:N
            switch data(j,i)
                case 3
                    red = [red, i];
                case 2
                    yellow = [yellow,i];
            end
        end
        
        % Branch plot
        sol = soln_struc.sol;
        sol = sol(any(sol,2),:); % truncates all-0 rows
        contparam = soln_struc.p(:,soln_struc.icp);
        contparam = contparam(any(sol,2),:); % truncates param values corresponding to all-0 rows in sol
        subplot('Position', [0.55, 0.55, 0.4, 0.35]);
        plot(contparam, vecnorm(sol'), 'LineWidth', 2);
        hold on;
        plot(contparam(1), norm(sol(1,:)), 'Marker', 'o', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
        plot(contparam(end), norm(sol(end,:)), 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
        plot(patterns(j,N+1), norm(patterns(j,1:N)), 'Marker', 'o', 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 10);
        hold off;
        xlabel('\mu'); ylabel('||u||');
        title('Bifurcation Diagram');
        
        % Graph plot
        subplot('Position', [0.05, 0.55, 0.4, 0.35]);
        h = plot(G, 'MarkerSize', 10, 'NodeColor', 'k', 'EdgeColor', 'k');
        highlight(h, red, 'NodeColor', [0.929 0.110 0.141]);
        highlight(h, yellow, 'NodeColor', [1.000 0.780 0.173]);
        title('Graph');
   
        % Grid plot
        subplot('Position', [0.07, 0.1, 0.9, 0.35])
        image(data);
        map = [0.000 0.000 0.000;
               1.000 0.780 0.173;
               0.929 0.110 0.141];
        colormap(map);
        colorbar;
        set(gca, 'YDir', 'normal');
        % labels and lines
        xticks(1:N)
        yticks(1:M)
        xl = arrayfun(@(x)xline(x,'w-','LineWidth',1),1.5:1:N-0.5);
        yl = arrayfun(@(y)yline(y,'w-','LineWidth',1),1.5:1:M-0.5);
        axis equal
        xlabel('Node'); ylabel('Pattern #');
        
        key = input("Enter 1 to go forwards, 0 to go backwards, or 'q' to quit: ");
    end
    close
end