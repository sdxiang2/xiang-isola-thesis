function [percentage_correct, percentage_unpredicted, cs_ipg, isp_ipg, cs_fpd, isp_fpd, degdist] = ConStats_analysis2(cs_output, isp_output)
% Compares conjectured isola predictions to continuation results.
%
% (in all of the following, N is the number of nodes and M is the number of
% graphs.)
% Input:
%   cs_output: NxMx2 array. (:,:,1) is the isola conjecture results.
%       (:,:,2) is the node degree.
%   isp_output: NxMx2 array. (:,:,1) is the continued isola data.
%       (:,:,2) is the node degree.
%
% Output:
%   percentage_correct: float. Percentage of SDNs conjectured to be isolas
%       that are actually isolas when continued.
%   percentage_unpredicted: float. Percentage of SDNs that are isolas that
%       were not predicted by the conjecture
%   cs_ipg: 1xM array. Number of conjectured isolas per graph.
%   isp_ipg: 1xM array. Number of continued isolas per graph.
%   cs_fpd: (N-1)x1 array. Fraction of SDNs that were conjectured to be
%       isolas, per graph.
%   isp_fpd: (N-1)x1 array. Fraction of SDNs that are isolas, per graph.
%   degdist: (N-1)x1 array. Empirical degree distribution.


% ---------- BEGIN CODE ----------

    cs_isolas = squeeze(cs_output(:,:,1));
    cs_degrees = squeeze(cs_output(:,:,2));

    isp_isolas = squeeze(isp_output(:,:,1));
    isp_degrees = squeeze(isp_output(:,:,2));

    N = size(cs_isolas,1);
    M = size(cs_isolas,2);

    % Mean degree for cs and isp
    % These should be identical
    cs_mean_degree = sum(cs_degrees,'all')/numel(cs_degrees);
    isp_mean_degree = sum(isp_degrees,'all')/numel(isp_degrees);
    disp([cs_mean_degree,isp_mean_degree]);

    % Percentage of conjectured isolas that are actually isolas
    isp_filter_cs_isola = isp_isolas(cs_isolas == 1);
    percentage_correct = sum(isp_filter_cs_isola,'all')/numel(isp_filter_cs_isola);
    fprintf('Percentage correct: %f\n', percentage_correct);

    % Percentage of actual isolas not predicted by the conjecture
    isp_cs_diff = isp_isolas - cs_isolas;
    percentage_unpredicted = sum(isp_cs_diff == 1,'all')/sum(isp_isolas == 1,'all');
    fprintf('Percentage unpredicted: %f\n', percentage_unpredicted);

    % Conjectured isolas per graph
    cs_ipg = sum(cs_isolas,1);

    % Actual isolas per graph
    isp_ipg = sum(isp_isolas,1);

    % Plotting isolas per graph
    figure
    histogram(cs_ipg, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'r');
    hold on
    histogram(isp_ipg, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'b');
    xlabel('# isolas');
    ylabel('# graphs');
    legend('Conjectured', 'Actual');
    hold off

    % Conjectured fraction isolas per degree
    maxdegree = max(cs_degrees,[],'all');
    cs_fpd = NaN(maxdegree+1,1);
    degdist = NaN(maxdegree+1,1);
    for i = 0:maxdegree
        filter_isolas = cs_isolas(cs_degrees == i);
        total = length(filter_isolas);
        num_isolas = nnz(filter_isolas(:) == 1);
        cs_fpd(i+1,1) = num_isolas/total;
        degdist(i+1) = total;
    end

    % Actual fraction isolas per degree
    isp_fpd = NaN(maxdegree+1,1);
    for i = 0:maxdegree
        filter_isolas = isp_isolas(isp_degrees == i);
        total = length(filter_isolas);
        num_isolas = nnz(filter_isolas(:) == 1);
        isp_fpd(i+1,1) = num_isolas/total;
    end

    % Plotting fraction isolas per degree   
    figure
    yyaxis left
    scatter(0:maxdegree, cs_fpd, 'o', 'filled', 'SizeData', 100);
    hold on
    scatter(0:maxdegree, isp_fpd, 'o', 'SizeData', 100);
    xlabel('Degree');
    ylabel('Fraction of isolas');
    yyaxis right
    scatter(0:maxdegree, degdist, '+', 'SizeData', 100);
    ylabel('Number of nodes');
    legend('Conjectured', 'Actual', 'Degree distribution', 'Location', 'northwest');

end