% generating random data
for i = 1:10
	n = 300;
	m = 300;
	C = rand(n,m);
	a = rand(n,1);
	b = rand(m,1);
	b = b / sum(b) * sum(a);%ensuring sum(a) = sum(b)
	save(sprintf("..\\datasets/randomC%d",i),'C')
	save(sprintf("..\\datasets/randoma%d",i),'a')
	save(sprintf("..\\datasets/randomb%d",i),'b')
end