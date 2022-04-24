function [comparisonArray,totals,byDegree] = FS_analysis(isp_output, isp_fsa, cs_fsa)
% Compares ISP and ConStats first-switch data.
%
% (in all of the following, N is the number of nodes and M is the number of
% graphs.)
% Input:
%   isp_output - NxMx2 array. (:,:,1) is the continued isola data (unused),
%       and (:,:,2) is the degree of each node.
%   isp_fsa - NxM cell array. Saved as "firstSwitchedArray" by
%       IsolaStatsParallel.
%   cs_fsa - NxM cell array. Saved as "firstSwitchedArray" by ConjectureStats.
%
% Output:
%   comparisonArray - NxM array. The number in the ith row and jth column
%       indicates the conjecture accuracy for node i in graph j. 3 means the
%       actual outcome matches the conjectured outcome. 2 means the actual
%       outcome(s) is a subset of the conjectured outcomes. 1 means neither of
%       these is true. -1 indicates an error.
%   totals - 1x4 array. The first number is the number of times the
%       conjecture was correct. The second is the number of times the
%       conjecture was partially correct (the actual first-switched was in
%       the set predicted by the conjecture). The third is the number of
%       times the conjecture was wrong. The fourth is the number of times
%       the analysis errored in some way.
%   byDegree - (maxdegree+1)x7 array. totals data broken down by degree. 
%       First column is degree. Second column is the number of nodes with 
%       that degree. Third column is the average number of conjectured
%       first-switches for nodes with that degree.


% ---------- BEGIN CODE ----------

if size(isp_fsa) ~= size(cs_fsa)
    error('Size of input arrays are not equal.')
end

degrees = isp_output(:,:,2);

[N,M] = size(isp_fsa);
comparisonArray = ones(N,M);

%gcp;

for i = 1:(N*M)
    isp = isp_fsa{i};
    cs = cs_fsa{i};
    if ~isnumeric(isp) || ~isnumeric(cs)
        comparisonArray(i) = -1;
    elseif isequal(isp,cs)
        comparisonArray(i) = 3;
    elseif all(ismember(isp,cs))
        comparisonArray(i) = 2;
    end
end

correct = sum(comparisonArray == 3, 'all');
subset = sum(comparisonArray == 2, 'all'); % isp subset of cs
incorrect = sum(comparisonArray == 1, 'all');
errored = sum(comparisonArray == -1, 'all');

fprintf('Correct: %g    Subset: %g  Incorrect: %g   Errored: %g', correct, subset, incorrect, errored);

totals = [correct subset incorrect errored]/(N*M);

maxdegree = max(degrees,[],'all');
byDegree = NaN(maxdegree+1, 7);
byDegree(:,1) = 0:maxdegree;
for d = 0:maxdegree
    filter_comparisonArray = comparisonArray(degrees == d);
    filter_cs_fsa = cs_fsa(degrees == d); % can you do this to cell arrays?
    mean_num_predictions = mean(cellfun('length',filter_cs_fsa));
    numdeg = length(filter_comparisonArray);
    byDegree(d+1,2) = numdeg; % number of nodes with degree d
    byDegree(d+1,3) = mean_num_predictions;
    byDegree(d+1,4) = length(filter_comparisonArray(filter_comparisonArray == 3))/numdeg; % correct
    byDegree(d+1,5) = length(filter_comparisonArray(filter_comparisonArray == 2))/numdeg; % subset
    byDegree(d+1,6) = length(filter_comparisonArray(filter_comparisonArray == 1))/numdeg; % incorrect
    byDegree(d+1,7) = length(filter_comparisonArray(filter_comparisonArray == -1))/numdeg; % error
end


%delete(gcp('nocreate'))

end