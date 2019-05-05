N = 300;
for n = 1:10
	noise = randn(2*N,1)*0.1;
	noise = reshape(noise, [N, 2]);
	primal = rand(N,1) * 2 * pi;
	primalpoints = [cos(primal) sin(primal)] + noise;
	source = primalpoints.*[1.3 0.9];
	target = primalpoints.*[0.9 1.1];
	C = zeros(N);
	for i = 1:N
		for j = 1:N
			C(i, j) = norm(source(i,:) - target(j,:))^2;
		end
	end
	save(sprintf('..\\datasets\\ellipsedataset%d.mat',n),'C') 
end