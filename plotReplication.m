function plotReplication(ts,hist,params)

l = params.l;
tlims = params.tlims;
epochs = params.epochs;

%% Plotting
% plot of growth trends
figure()
surf(l,ts,hist,'LineStyle','none')
colormap hot
colorbar
hold on
%rskip = round(linspace(1,size(hist,1),25));
%cskip = round(linspace(1,size(hist,2),25));
%surf(l(cskip),ts(rskip),hist(rskip,cskip),'FaceColor','none','EdgeColor','white');
hold off

title('Growth Trend of Genomes, seperated by Length')
zlabel('Number of Genomes')
ylabel('Time')
xlabel('Length')
%legend(strread(num2str(l),'%s'), 'Location', 'BestOutside')
%%
% plot of average length
times = ismember(ts,[0,(1:epochs)*tlims(2)]);
averages = (l*hist')'./sum(hist,2);
figure()
plot(ts,averages,'-',ts(times),averages(times),'o')
title('Average Genome Length  over Time')
xlabel('Time')
ylabel('Length')
%%
% plot of number of genomes (sanity check)
figure()
plot(ts,sum(hist,2))
title('Number of Genomes over Time')
xlabel('Time')
ylabel('Number')
%%
% plot of t_final distributions
figure()
plot(l,hist(end,:),'o')
title('Distribution of Lengths at t_{final}')
xlabel('Length')
ylabel('Number of Genomes')