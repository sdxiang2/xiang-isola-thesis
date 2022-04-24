function ConjectureStats2(p,generator,savestruc)
% Generates a series of random graphs, and for each node/SDN, determines
% whether it should result in a figure 8 isola based on our conjecture. The
% conjecture is that an SDN will result in an isola if k = argmin(cj) is
% unique and d_k < d_SDN.
% 
% Inputs:
%   p - [0,1] (probability of edge for ER graphs) or integer (initial graph
%       size/number of edges added per node for BA graphs)
%   generator - {'ER','BA'}. Specifies random graph generator
%   save_struc - boolean. If true, saves the graph laplacians of the
%       generated graphs.
%
% Outputs:
%   None. Saves isola-prediction data and graph laplacians.


% ---------- BEGIN CODE ----------

    gcp;
    tic
    M = 200; % # of graphs
    N = 50; % # of nodes per graph

    %M = 3;
    %N = 10;
    
    if strcmp(generator,'ER')
        GL = ER_GL_generator(M,N,p,0);
    elseif strcmp(generator,'BA')
        GL = BA_GL_generator(M,N,p,0);
    else
        error("Error: invalid generator. Should be either 'ER' or 'BA'");
    end

    if savestruc
        %filename = sprintf('/users/sxiang2/data/sxiang2/ConStats_data/GL_%g.mat',p);
        filename = sprintf('/Users/stacey/Documents/MATLAB/Oscar_files/graphlaplacians/testConStatsGL_%g.mat',p);
        save(filename, 'GL');
    end

    % For SDN N in graph M:
    isIsolaArray = NaN(N,M); % do we think it'll be an isola?
    degreeArray = NaN(N,M); % what is its degree?
    firstSwitchedArray = cell(N,M); % what do we think will be the first node that switches?

    parfor i = 1:(N*M)
    %for i = 1:(N*M)
        % Getting graph number and SDN number from i
        graph_num = ceil(i/N);
        node_num = mod(i-1,N) + 1;
        %disp([graph_num,node_num]);
        % Getting graph laplacian
        L = squeeze(GL(graph_num,:,:)); 
        neighbors = find(L(node_num,:));
        neighbors(neighbors == node_num) = [];
        degreeArray(i) = length(neighbors);
        min_cj = inf;
        argmin_cj = [];
        for n = neighbors
            degree = -L(n,n);
            second_neighbors = find(L(n,:));
            second_neighbors(second_neighbors == n) = [];
            shared_neighbors = intersect(neighbors,second_neighbors);
            cj = degree - length(shared_neighbors); % the node with min cj should switch first
            if cj < min_cj
                argmin_cj = n;
                min_cj = cj;
            elseif cj == min_cj
                argmin_cj = [argmin_cj,n];
            end
        end

        if length(argmin_cj) == 1 && -L(argmin_cj,argmin_cj) < -L(node_num,node_num)
            isIsolaArray(i) = 1;
        else
            isIsolaArray(i) = 0;
        end

        firstSwitchedArray{i} = argmin_cj;
    end

    output(:,:,1) = isIsolaArray;
    output(:,:,2) = degreeArray;
    
    %filename = sprintf('/users/sxiang2/data/sxiang2/ConStats_data/ConStats_%gg_%gn_%g.mat',M,N,p);
    filename = sprintf('/Users/stacey/Documents/MATLAB/Oscar_files/constatsdata/ConStats_output_%g.mat',p);
    save(filename, 'output');
    filename = sprintf('/Users/stacey/Documents/MATLAB/Oscar_files/constatsdata/ConStats_fsa_%g.mat',p);
    save(filename, 'firstSwitchedArray');
    
    delete(gcp('nocreate'))
    toc
end












