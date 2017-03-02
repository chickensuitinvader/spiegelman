function theOtherPlottingFunction(ts,hists,params,gamma)
l = params.l;
times = [0,(1:params.epochs)*params.tlims(2)];
for q = 1:length(gamma)
    gString = sprintf('Transfer Proportion = %.3f',gamma(q));
    for p = 1:length(params.kappa)
        a = find(ismember(ts{p,q},times));
        a = a(mod(1:length(a),2) == 1);
        t(:,p,q) = ts{p,q}(a);
        average = (l*hists{p,q}')'./sum(hists{p,q},2);
        y(:,p,q) = average(a);
    end
    
    % plots of average lengths over time
%     figure()
%     plot(t(:,:,q),y(:,:,q));
%     h = legend(num2str(params.kappa(:)),'Location','BestOutside');
%     v = get(h,'title');
%     set(v,'string','Transfer Parameter Value');
%     title({'Mean Genome Length After Transfer Step',gString});
%     xlabel('Time')
%     ylabel('Length')
    
end

temp(:,:) = y(end,:,:);
figure()
plot(params.kappa,temp)
title('Mean Genome Length at t_{final}')
xlabel('Transfer Parameter Value')
ylabel('Length')
h = legend(num2str(gamma(:)),'Location','BestOutside');
v = get(h,'title');
set(v,'string','Transfer Percentage');

figure()
surf(gamma,params.kappa,temp,'EdgeColor','interp')
title('Mean Genome Length at t_{final}')
xlabel('Transfer Proportion')
ylabel('Transfer Parameter Value')
zlabel('Length')
end

