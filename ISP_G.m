function [new_isp_fsa] = ISP_G(comparisonArray,GL,isp_fsa,ds)
% Finds SDNs for which the predicted first-switch is incorrect, runs them
% again with smaller step size, returns an updated isp_fsa for use in
% FS_analysis
%
% Inputs:
%   comparisonArray - NxM array. Output by FS_analysis.
%   GL - MxNxN array. Graph laplacians. Saved by ConjectureStats.
%   isp_fsa - NxM cell array. Saved by IsolaStatsParallel. Entry at (n,m)
%       are the actual first switches for SDN n, graph m when continued.
%   ds - float. New step size. Should be smaller than 0.001 (original step
%       size).
%
% Outputs:
%   new_isp_fsa - NxM cell array. isp_fsa with incorrectly predicted first
%       switches replaced with new result when run with stepsize ds.


% ---------- BEGIN CODE ----------
    
    gcp;
    tic

    M = size(GL,1); % # of graphs
    N = size(GL,2); % # of nodes

    filter_ones = comparisonArray == 1;
    incorrects = find(filter_ones);

    % Continuation options
    pars = [0.01,0.4];
    u_ = sqrt(1 - sqrt(1 - pars(2)));
    options = [0,0.5];

    % Continuation
    temp = cell(size(incorrects));
    parfor i = 1:length(incorrects)
    %for i = 1:length(incorrects)
        graph = ceil(incorrects(i)/N);
        node = mod(incorrects(i)-1,N) + 1;
        L = squeeze(GL(graph,:,:));
        sh = @(u,p) LDS_RHS(u,p,L);
        u0 = zeros(N,1);
        u0(node) = u_;
        [ss,os] = TangentContinuationLDS(L,sh,u0,pars,2,ds,(0.001/ds)*10000,options);
        pat = GetPatterns(ss,os);
        temp{i} = FindFirstSwitched(pat);
    end

    new_isp_fsa = isp_fsa;
    new_isp_fsa(incorrects) = temp;
    
    delete(gcp('nocreate'))
    toc

end