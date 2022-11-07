Files = dir('MPC_*');
idx = 1601;
XData = 0:0.2:320;
YData_MPC = [];
horizon = [];
for i = 1:length(Files)
    horizon(i) = str2num(Files(i).name(13:14));
    fig = open(Files(i).name);
    YData_MPC(i,:) = fig.Children(2).Children(2).YData(1:idx);
    if i == 1
        YData_OCP = fig.Children(2).Children(1).YData(1:idx);
    end
    close(fig)
end

plot(horizon, mean(abs(YData_MPC - YData_OCP),2), LineWidth=2)
xlabel("Prediction horizon in $s$")
ylabel("Mean absolute error")