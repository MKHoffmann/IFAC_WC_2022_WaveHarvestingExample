fig = figure();
% out = arrayfun(@(i) [i.XData, i.YData], ax.Children, 'uni', 0);
dat = vertcat(out{:});
dat_s = sortrows(dat);
plot(dat_s(2:end-1,1), dat_s(2:end-1,2), 'x', LineWidth=2, MarkerSize=8)
hold on
plot(dat_s([1 end],1), dat_s([1 end],2), 'o', LineWidth=2, MarkerSize=8)
plot(-16.2891, 0.465273, '*', LineWidth=2, MarkerSize=8)
hold off
ylim([-0.1, 5])
grid on
ax = gca;
ax.LineWidth = 1;
ax.GridAlpha = 1;
ax.GridColor = 0.65*ones(1,3);
legend("valid weights", "extreme points", "weight controller", fontsize=12)
ax.FontSize = 13;
ax.YAxis.TickValues = 0:0.5:5;
xlabel("Extracted energy in MJ")
ylabel("Damage")
saveas(fig, "mpc_paretofront_9.svg")