function d = calcDeltaM( l, ps, mode )
%calcDeltaM Calculates the probability that a combination of point mutations
%will map from i to j, ignores any mappings for i to i
%   l = array of lengths
%   ps = array of probailities
%   d = matrix of point mutation probabilities

if ~exist('mode','var'); mode = 1; end;

d = eye(length(l));

if mode == 0 % No mutation
    return
end

for i = 1:length(l) %FROM
    dif = l-l(i);
    d(i,:) = (bindif(dif,l(i),ps))';
end
end

