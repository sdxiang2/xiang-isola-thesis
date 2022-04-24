function[Summary_Stats] = erdos_renyi_distribution(p)
%% Setup
n = 50;
%p = 0.4;

M = 100000;
D = [1:1:n-1]; % S: degrees

Summary_Stats = zeros(length(D),4);
Summary_Stats(:,1) = D';
fprintf('degree IsolaGuess IsolaGuessPred\n')
ell = 1;

%% MCMC simulations
for d = D
	P2 = 0; % counter for conjectured isolas

	%	a = randsample(n-1,d)'+1;   % bottleneck of computation
	a = [1:d]+1; % S: nodes 2 to d+1. We can do this because it doesn't matter which d nodes we connect to.

	vals = [];
	for j=1:M
        % S: creating random Erdos-Renyi adjacency matrix
		G = rand(n,n) < p;
		G(1,:) = 0;
		G(1,a) = 1;
		G = triu(G,1);
		G = G+G';
		s = sum(G); % degrees

		% P2=P(deg(k)<d: k=argmin(bar{c}_j: j is neighbor of SDN))
		%		[~,k] = min(sum(G(a,G(1,:)==0),2));
		[~,k] = min(sum(G(a,[d+2:end]),2));
        % S: This does the min(X1,...,Xd) thing from the probability
        % distribution notes
		if s(a(k))<d % if the degree of the argmin node is less than the degree of the SDN
			P2 = P2+1;
		end

		Ideg = sum(G(a,[d+2:end]),2); % |N(k) intersect I|
		[m,k] = min(Ideg);
		K = a(Ideg==m); % all nodes with minimum Ideg
		vals = [vals; length(K)]; % number of nodes with minimum Ideg	
	end

	Fd = cdfdeg([0:d-1],d,p); % S: cdf of connecting to [0:d-1] nodes out of d given p
	Gd = 1-(1-cdfdeg([0:n-d-1],n-d,p)).^d; % S: cdf of complement of above
	Pd = [Fd(1),diff(Fd)]; % S: pdf of connecting to [0:d-1] nodes out of d given p
	Qd = [Gd(1),diff(Gd)]; % S: pdf of complement of above
	Z = conv(Pd,Qd); % S: not sure I understand why we have to convolute these two?
	PP2 = sum(Z(1:d-1)); % S: cdf of Z

	Summary_Stats(ell,2:end) = [P2/M, PP2, mean(vals)];
	fprintf('%d  %.4f  %.4f\n', d, P2/M, PP2)
	ell = ell+1;
end

%% Plot results
figure(1)
subplot(1,2,1)
plot(Summary_Stats(:,1)/n,Summary_Stats(:,2),'r*', Summary_Stats(:,1)/n,Summary_Stats(:,3),'b-o');
grid
subplot(1,2,2)
plot(Summary_Stats(:,1)/n,Summary_Stats(:,end),'r*');
grid


%% Cumulative distribution of degree d on graph with n vertices
function y = cdfdeg(d,n,p)
	y = cdf('Binomial',d,n-1,p);
end

end
