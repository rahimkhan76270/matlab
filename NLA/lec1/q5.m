n=8;

coeffs=randperm(n); % generate random values
disp(coeffs);
disp(length(coeffs));% calculate length of coeffs
disp(8-1:-1:0); % length-1 to 0 all values 
disp((length(coeffs)-1:-1:0).*coeffs); % take elementwise product
