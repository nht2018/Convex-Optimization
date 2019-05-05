for i = 0:9
	for j = 0:9
		temp = xlsread(sprintf('images/image%d_%d.csv',[i,j]));
		ker = [1 1 0;1 1 0;0 0 0];
		temp = imfilter(temp,ker,'conv');
		result = temp(2:2:32,2:2:32);
        result = reshape(result,[256,1]);
        result = result / sum(result) * 128*16^2;
		save(sprintf('..\\datasets/image%d_%d',[i,j]),'result')
	end
end
