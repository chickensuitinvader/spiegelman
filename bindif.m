function [ v ] = bindif(z,n,p)
% returns a row vector corresponding to the probabilites of changes, z, of
% something with size n, with probability of positive change pa, and
% negative change pd
% If z is discretised (non-unit intervals), then returns a cumulative sum
% of pdfs around the discretisations

% assumption that p is in form [pa,pd,pba,pbd]
z0 = [2*z(1)-z(2),z,2*z(end)-z(end-1)];
z1 = repmat([-Inf,z0(1):z0(end),Inf],n+1,1);
p1 = p(2-(z1>0));
p2 = p(1+(z1>0));
x = repmat((0:n)',1,size(z1,2));
pdf = sum(binopdf(x+abs(z1),n,p1).*binopdf(x,n,p2));
zero = pdf(z1(1,:) == 0);
cdf = cumsum(pdf);
v = diff(cdf(ismember(z1(1,:),z0)));
if sum(z < 0) >= 1
    ind = find(z==0);
    er = v(ind) - zero;
    v(ind) = v(ind) - er;
    v(ind-1) = v(ind-1) + er;
end
v = v(1:end-1);


