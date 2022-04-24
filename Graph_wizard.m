function [neighbors] = Graph_wizard(L,node,display)
    if display
        fprintf('This node has degree %d.\n', -L(node,node));
    end
    neighbors = find(L(node,:));
    neighbors(neighbors == node) = [];
    if display
    fprintf('The neighbors of this node are:\n');
        for n = neighbors
            degree = -L(n,n);
            second_neighbors = Graph_wizard(L,n,false);
            shared_neighbors = intersect(neighbors,second_neighbors);
            value = degree - length(shared_neighbors); % the node with minimal value should switch first
            fprintf('%d \t degree: %d \t shared neighbors: %d \t value %d\n',...
                n, degree, length(shared_neighbors), value);
        end
    end
end