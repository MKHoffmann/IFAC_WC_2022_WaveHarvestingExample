% 1. Run the OPC with large time Horizon
%       - Find global weight for weighted sum.
%       -
global PathToParameters
PathToParameters= 'C:\Users\MKH\Lennart\StochasticWavesMOOP\src\PolySurge_inputs.mat';
load(PathToParameters);
filenameMOOP = ['MOOPStochastic_400seconds.mat'];
%%
timehorizon = 400;                          % shoud be self explanatory
timestep = 0.2;                             % shoud be self explanatory
dmg_at_horizon = 8;
nHorizon = round(timehorizon/timestep);
% Number of Pareto Points
time = [0:timestep:(nHorizon)*timestep];    % Create array with discrete time steps
ocp_t = time;


[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(timehorizon,timestep);
monoW = monochromaticWave();
dval = arrayfun(@(t) FBMStochasticWave(t),[0:timestep:((d.length()-1)*timestep)]);
ocp.set_value(d,dval);

costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);
ocp.minimize( costfun*p_params );
% [p_params, ep_ocp, w_ep] = scalarize_moocp( ocp, costfun, method='ws', normalize='fix' );

tic
nPoints = 15;
w1= linspace(0.05, 1-0.05,nPoints);
w = [w1 ;1-w1]';

Point = ceil(nPoints/2+2);
sol = [];
ResultsMOOP = cell(nPoints,1);
for i = 1:nPoints
    ocp.set_value(p_params,w(i,:));
    if ~isempty(sol)
        if i == 2
            solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
                'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
                'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
            ocp.solver('ipopt', struct(), solver_opt)
        end
%         ocp.set_initial(x(:,2:end),sol(i-1).value(x(:,2:end)))
        ocp.set_initial([ocp.x; ocp.lam_g], sol(i-1).value([ocp.x; ocp.lam_g]))
    end
    sol = [sol ocp.solve()];
    MOOPSolution = struct;
    MOOPSolution.weigths = w(i,:);
    MOOPSolution.x = ocp.value(x);
    MOOPSolution.u=ocp.value(u);
    MOOPSolution.ParetoPoint=[ocp.value(x(6,end)) ocp.value(x(7,end))];
    ResultsMOOP{i} = MOOPSolution;
end
%% Chose Point here! Important choice. 
Point = 6;
% Validate with better integrator. 
d = sol(Point).value(d);


U = @(t) interp1(time,sol(Point).value(u),t,'linear');  %% Zero order hold
validationdgl = @(t,x) [Ac * x(1:5) - Bc * 1e6 * U(t) * gamma * x(2) + Bc * FBMStochasticWave(t)];
k = ode89(validationdgl,[0,timehorizon],x0(1:5));
f2 = figure(2)
plot(k.x,k.y(2,:))
hold on 


plot(time,sol(Point).value(x(2,:)))

RMSE = SignalDifference(k.x,k.y(2,:),time,sol(Point).value(x(2,:)))


%% Invesigate the damage over the horizond
f3 = figure(3)
for i = 1:nPoints
hold on 
plot(ocp_t,ResultsMOOP{i}.x(7,:))
end

%% Chose damage at the end of horizond 
try
    clear ResultsMPC;
end
DamageTarget = 0.35;
timehorizon = 3000;
DamageAtIncrement = @(t) interp1([0 timehorizon],[0 DamageTarget],t,"linear");


Point = 6; %StartingPoint
% After confirming the time step is adequate start MPC here

MPCtimehorizon = 60;
Seed = 14;
nPoints = 15;
timestep = 0.5;
AdjustPointInterval = 50;


filename = ['TrackingDamageIncremental_Faster_Switch_ON_Down_DamageTarget_' num2str(DamageTarget) '_Seed_' num2str(Seed) '.mat'] 
w1= linspace(0.05, 1-0.05,nPoints);
w = [w1 ;1-w1]';
% filename = ['Stochastic,' datestr(now,'DD_HH_MM') '_MPC_Horizon_weightControl' num2str(MPCtimehorizon) 'Point_7.mat'];
x0 =   [0.0515     0.1392   314.4363   94.1062  190.4844         0         0   19.6030   39.2764  -49.7962       0         0]';

[ocp,x,u,d,x0,x0_p] = initializeOCPENERGY(MPCtimehorizon,timestep,x0=x0);
costfun = ([x(6,end) x(7,end)]);
p_params = ocp.parameter(2,1);
ocp.set_value(p_params, w(Point,:))
ocp.minimize( costfun*p_params );
starttime = [0:0.5:timehorizon];
AppliedSignal=[];
cost_hist = [];
PointHist = [];
last_i = 0;
f1 = figure(144)


plot([0:5:timehorizon],arrayfun(@(t) DamageAtIncrement(t),[0:5:timehorizon]))
title(filename)
colormap(f1,"cool")
 cb = colorbar; 
 cb.Label.String = 'Pareto Point';

caxis([1,nPoints]); 
hold on 
for i = 1:length(starttime)
    nHorizon = round(MPCtimehorizon/timestep);     % Number of Pareto Points
    time = [starttime(i):timestep:(nHorizon)*timestep+starttime(i)];    % Create array with discrete time steps

    if i == 2
        solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-8, 'warm_start_mult_bound_push',1e-10, ...
            'warm_start_slack_bound_push', 1e-10, 'warm_start_bound_push', 1e-8, ...
            'warm_start_bound_frac',1e-8, 'warm_start_slack_bound_frac',1e-10, ...
            'print_level', 1);
        ocp.solver('ipopt', struct(), solver_opt)
    end
    ocp.set_value(d,arrayfun(@(t) FBMStochasticWave(t,"Seed",Seed),time));

    % if not the first MPC step
    if ~(i == 1)
        % reset solver for warmstarting
        ocp.set_value(x(:,1), solOld.value(x(:,2)));
%         ocp.set_value(x(6:7,1), zeros(2,1));
        ocp.set_value(u(1),solOld.value(u(2)));

        ocp.set_initial([ocp.x; ocp.lam_g], solOld.value([ocp.x; ocp.lam_g]))
        %         ocp.set_initial(x(:,2:end-1),solOld.value(x(:,1:end-2)));
        %         ocp.set_initial(u(2:end-1),solOld.value(u(1:end-2)));
    end

    solOld = ocp.solve();
    solMpc = struct;
    solMpc.x = solOld.value(x);
    solMpc.u = solOld.value(u);
    solMpc.time = solOld.value(time);
    cost_hist(:,i) = solOld.value(x(6:7,2));
    PointHist(i) = Point;
    if ((mod(i, AdjustPointInterval) == 0))
        dmg_per_s = (cost_hist(2,end)-cost_hist(2,end-(AdjustPointInterval-1)))/(AdjustPointInterval*timestep);
        dt = timehorizon - starttime(i)+timestep;
        current_dmg = (cost_hist(2,end));
        if sum((AppliedSignal(2,end-(AdjustPointInterval-2):end))) < 50
            Point = Point; % Too little activity to justify point change 
        elseif current_dmg + dmg_per_s*dt > DamageTarget && Point(end) > 1 %% und threshholdspannung
            Point = Point - 1;
            ocp.set_value(p_params,w(Point,:));
        elseif current_dmg + dmg_per_s*dt < 0.8*DamageTarget && Point(end) < 15 && ((mod(i, AdjustPointInterval*2) == 0))
            Point = Point + 1;
            ocp.set_value(p_params,w(Point,:));


        end

    scatter(solOld.value(time(1)),solOld.value(x(7,1)),120,Point,'filled') 
    disp('NewPoint = ')
    disp(num2str(Point))
    end
    AppliedSignal = [AppliedSignal [time(1); solOld.value(u(1))]];
    ResultsMPC(i) = solMpc;

end

%
% f4 = figure(4)
% for i = 1:length(starttime)
%
%
% Divergence = ResultsMPC{i}.u-arrayfun(@(t) OCP_U(t),ResultsMPC{i}.time);
% plot(ResultsMPC{i}.time,Divergence)
% hold on

% f5 = figure(5)
% plot(AppliedSignal(2,:)-arrayfun(@(t) OCP_U(t),[0:timestep:AppliedSignal(1,end)]))


 save (filename,'ResultsMPC', 'PointHist', 'DamageAtIncrement', 'DamageTarget','cost_hist')
 %% TRY out different stuff to get turnpike visualized
% i = 10;
% f5 = figure(12)
%
% plot(ResultsMPC{i}.time,(ResultsMPC{i}.u-arrayfun(@(t) OCP_U(t), ResultsMPC{i}.time)).^2)
% hold on
%









function [RMSE]  = SignalDifference(t1,sig1,t2,sig2)
sig2 = @(t) interp1(t2,sig2,t);
RMSE = 0;
for i = 1:length(t1)
    RMSE = RMSE + (sig1(i)-sig2(t1(i)))^2;



end
RMSE = sqrt(RMSE/length(t1));

if RMSE >= 0.02
    warning('Please reconsider the time step! might be to large')
end
end