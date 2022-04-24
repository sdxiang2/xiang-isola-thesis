 function [isIsola] = CheckIsola(patterns)
% From a bifurcation diagram's pattern data, determines whether the
% bifurcation diagram is a figure 8 isola.
%
% Inputs:
%   patterns - Pattern data generated from GetPatterns.
% Outputs:
%   isIsola - boolean. Whether or not the bif. diagram is an isola.


% ---------- BEGIN CODE ----------

    isIsola = false;

    [~,SDN] = max(patterns(1,:));

    n = size(patterns,2) - 1;
    
    if size(patterns,1) < 4
        return
    end

    num_uminus = 0;
    for r = 1:n
        if abs(patterns(1,SDN) - patterns(4,r)) < 0.1
            num_uminus = num_uminus+1;
        end
    end
    
    if size(patterns,1) > 3 && abs(patterns(1,SDN) - patterns(4,SDN)) < 0.1 && num_uminus == 2
        isIsola = true;
    end
    
 end