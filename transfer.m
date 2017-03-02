function n = transfer(l,n,g,k)
ws = l.^k.*n;
n = g*sum(n)*(ws/sum(ws));
end