function [GL] = BA_GL_generator(numgraphs,numnodes,k,seed)
    rng(seed);
    GL = zeros(numgraphs,numnodes,numnodes);
    connected_graph = ones(k);
    connected_graph = connected_graph - k*eye(k);
    initial_graph = zeros(numnodes);
    initial_graph(1:k,1:k) = connected_graph;
    for g = 1:numgraphs
        graph = initial_graph;
        for r = k+1:numnodes
            degrees = -diag(graph(1:r-1,1:r-1));
            new_edges = zeros(k,1);
            for e = 1:k
                p = degrees/sum(degrees);
                cdf = cumsum(p);
                U = rand;
                for i = 1:r-1
                    if U < cdf(i)
                        new_edges(e) = i;
                        degrees(i) = 0;
                        break
                    end
                end
            end
            graph(r,new_edges) = 1;
            graph(new_edges,r) = 1;
            for e = 1:k
                node = new_edges(e);
                graph(node,node) = graph(node,node) - 1;
            end
            graph(r,r) = -k;
        end
        GL(g,:,:) = graph;
    end
end