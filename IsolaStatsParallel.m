function IsolaStatsParallel(p,generator,save_struc)
% Generates a series of random graphs, and for each node/SDN,
% determines whether it results in a figure-8 isola. Note that this code
% won't run for Barabasi-Albert graphs if N/k is less than 2.
%
% Inputs:
%   p - [0,1] (probability of edge for ER graphs) or integer (initial graph
%       size/number of edges added per node for BA graphs)
%   generator - {'ER','BA'}. Specifies random graph generator
%   save_struc - boolean. If true, saves soln_struc, opt_struc, and pattern
%       data.
%
% Outputs:
%   None. Saves isola-check data and soln_struc, opt_struc, and pattern
%   data.


% ---------- BEGIN CODE ----------
    
    gcp;
    tic
    M = 200; % # of graphs
    N = 50; % # of nodes per graph

    %M = 10;
    %N = 15;
    
    % Generating graphs
    if strcmp(generator,'ER')
        GL = ER_GL_generator(M,N,p,0);
    elseif strcmp(generator,'BA')
        GL = BA_GL_generator(M,N,p,0);
    else
        error("Error: invalid generator. Should be either 'ER' or 'BA'");
    end
    
    % Continuation options
    pars = [0.01,0.4];
    u_ = sqrt(1 - sqrt(1 - pars(2)));
    options = [0,0.5];
    % For SDN N in graph M:
    isIsolaArray = NaN(N,M); % does it result in an isola?
    degreeArray = NaN(N,M); % what is its degree?
    firstSwitchedArray = cell(N,M); % what is the first node(s) that switches?
    
    % Continuation
    parfor i = 1:(N*M)
        graph = ceil(i/N);
        node = mod(i-1,N) + 1;
        L = squeeze(GL(graph,:,:));
        sh = @(u,p) LDS_RHS(u,p,L);
        u0 = zeros(N,1);
        u0(node) = u_;
        [ss,os] = TangentContinuationLDS(L,sh,u0,pars,2,0.001,10000,options);
        pat = GetPatterns(ss,os);
        if save_struc
            parsave(ss,os,pat,i,N);
        end
        isIsolaArray(i) = CheckIsola(pat);
        degreeArray(i) = -L(node,node);
        firstSwitchedArray{i} = FindFirstSwitched(pat);
    end
    
    output(:,:,1) = isIsolaArray;
    output(:,:,2) = degreeArray;
    
    filename = sprintf('/users/sxiang2/data/sxiang2/ISP_data/ISP_output_%g.mat',p);
    %filename = '/Users/stacey/Documents/MATLAB/Oscar_files/test_ISPvsCS/testISP_output.mat';
    save(filename, 'output');
    filename = sprintf('/users/sxiang2/data/sxiang2/ISP_data/ISP_fsa_%g.mat',p);
    %filename = '/Users/stacey/Documents/MATLAB/Oscar_files/test_ISPvsCS/testISP_fsa.mat';
    save(filename, 'firstSwitchedArray');
    
    delete(gcp('nocreate'))
    toc
end

function parsave(soln_struc,opt_struc,patterns,i,N)
    graph = ceil(i/N);
    node = mod(i-1,N) + 1;
    filename = sprintf('/users/sxiang2/scratch/ISP_data/graph%dSDN%d.mat',graph,node);
    %filename = sprintf('/Users/stacey/Documents/MATLAB/Oscar_files/test_ISPvsCS/test_parsave_data/graph%dSDN%d.mat',graph,node);
    save(filename,'soln_struc','opt_struc','patterns');
end

        
    