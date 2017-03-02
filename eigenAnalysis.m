kappas = linspace(-2,2,100);
close all
M = cell(length(kappas),1);
V = M;
S = M;
figure()
for i = 1:length(kappas)
   M{i} = b*d*diag((1./params.l).*(params.l.^kappas(i)));
   [V{i},S{i}] = eigs(M{i});
   plot(params.l,V{i}*S{i})
   title(num2str(kappas(i)));
   legend(num2str(diag(S{i})));
   ylim([-1 1]);
   pause(0.125);
end