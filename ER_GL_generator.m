function [GL] = ER_GL_generator(numgraphs,numnodes,p,seed)
 % Generates a set of Erdos-Renyi graph laplacians.
 %
 % Inputs:
 %  numgraphs - int. Number of graphs to generate.
 %  numnodes - int. Number of nodes in each graph.
 %  p - [0,1]. Probability that any given edge will appear in the graph.
 %  seed - Random seed used.
 % Outputs:
 %  GL - numgraphs x numnodes x numnodes matrix containing the set of graph
 %      laplacians.

    rng(seed);
    GL = zeros(numgraphs,numnodes,numnodes);
    for g = 1:numgraphs
        for r = 2:numnodes
            for c = 1:r-1
                edge = rand <= p;
                GL(g,r,c) = edge;
                GL(g,c,r) = edge;
            end
        end
        for r = 1:numnodes
            GL(g,r,r) = -sum(GL(g,r,:));
        end
    end
end