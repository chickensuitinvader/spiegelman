function dndt = replicationsystem(t,n,l,d,b,k,r,rj,MODE)

% reshape
n = n(:);
l = l(:);

% if MODE == 0
%     dndt = (n/sum(n)*k*r.*(l>rj))./l;
% elseif MODE == 1
%     dndt = (n/sum(n)*k*r.*(l>rj))./l - mu(d).*n./l + eta(d,l,n);
% else
%     dndt = (n/sum(n)*k*r.*(l>rj))./l - mu(d).*n./l + eta(d,l,n) - mu(b).*n./l + eta(b,l,n);
% end

m = b'*d';
dndt = m*((n/sum(n)*k*r)./l).*rj;

end