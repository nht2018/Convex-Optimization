N = 300;
for n = 1:10
	i = 1;
	primalpoints = zeros(N, 2);
	while i <= N
		temp = rand(1,2);
		if norm(temp) <= 1
			primalpoints(i,:) = temp;
			i = i+1;
		end
	end
	target = primalpoints;
	target(:,1) =  (abs(primalpoints(:, 1))+2).*sign(primalpoints(:, 1));
	source = primalpoints;
	C = zeros(N);
	for i = 1:N
		for j = 1:N
			C(i, j) = norm(source(i,:) - target(j,:))^2;
		end
	end
	save(sprintf('..\\datasets\\Caffarellidataset%d.mat',n),'C')
end