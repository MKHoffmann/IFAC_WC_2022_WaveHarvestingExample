num = 12;
t_bound = [82 120];
file_u = "MPC_Horizon_" + num + ".fig";
file_theta = "ThetaAndThetaDotMPC_OCP_" + num + "Seconds.fig";
fig_theta = open(file_theta);
xdata = fig_theta.Children(2).Children(1).XData;
is_in_t = xdata >= t_bound(1) & xdata <= t_bound(2);
theta_mpc = [fig_theta.Children(2).Children(3).XData(is_in_t); fig_theta.Children(2).Children(3).YData(is_in_t)];
theta_gt = [fig_theta.Children(2).Children(4).XData(is_in_t); fig_theta.Children(2).Children(4).YData(is_in_t)];
thetadot_mpc = [fig_theta.Children(2).Children(1).XData(is_in_t); fig_theta.Children(2).Children(1).YData(is_in_t)];
thetadot_gt = [fig_theta.Children(2).Children(2).XData(is_in_t); fig_theta.Children(2).Children(2).YData(is_in_t)];
close(fig_theta);

fig = open(file_u);
fig.Children(1).FontSize=12;
fig.Children(1).String{3} = 'Wave Excitation in Nm';
fig.Children(2).Children(1).LineWidth=2;
fig.Children(2).Children(2).LineWidth=2;
ylim([0, 680])
ax = gca;
ax.XAxis.FontSize=12;
ax.YAxis.FontSize=12;
ax.YAxis.Label.String = '$u$ in \SI{}{\kilo\volt^2}';
xlim(t_bound);
fig.Children(1).Location = 'northwest';
fig.Position = [934   585   954   607];
hold on
yyaxis right
ax.YAxis(2).FontSize=12;
ax.YAxis(2).Color = 'k';
ax.YAxis(2).Label.String = ["Angle in \SI{}{\degree}";"Angle velocity in \SI{}{\degree\per\second}"];
plot(theta_mpc(1,:), theta_mpc(2,:), 'LineWidth', 2);
plot(thetadot_mpc(1,:), thetadot_mpc(2,:), '-','LineWidth', 2, Color=[0.4660 0.6740 0.1880]);
ylim([-40, 40])
fig.Children(1).String{3} = '$theta$ for MPC';
fig.Children(1).String{4} = '$delta$ for MPC';
% yyaxis left
% hold off