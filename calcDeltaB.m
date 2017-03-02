function b = calcDeltaB(l,ps,mode)
if ~exist('mode','var'); mode = 2; end;
b = eye(length(l));
if abs(mode) >= 2
    for i = 1:length(l)
        b(i,:) = ((1-binopdf(0,l,ps(3))).*(1-binopdf(0,l,ps(4))))/(2i);
        big = (1:length(l)) > 2i;
        b(i,big) = 0;
    end
    b = b + diag(1-sum(b-diag(diag(b))));
end

