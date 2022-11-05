


%  Define MPC parameters. These can and should be chosen by the user.
% Suggested
% values

MPCtimehorizon  = 60;          % [1-70]  MCP Timehorizon. After ~70 seconds the solution does not improve.
timehorizon     = 400;         % [1-inf]  How long
timestep        = 0.5;         % [0.05-1] MPC timestep, i.e. the discretisation of the ocp.
%          A step larger then 1 is not recommended.
Seed            = 2  ;         % [1-10]   Seed of the Wave distrurbance. Seeds [1-10] have been provided.
Damagereduction = 0.4;         % [0-1]    This implementation uses Multicriterial Optimization.
%          Therefore wen need to specify a weight for
%          the damage contribution

plotting        = true;        % If plotting the solution of the MPC is shown as it progresses.
% but slows down the algorithm!


load('PolySurge_inputs.mat');
load ('WaveData.mat')
nSteps          = round(MPCtimehorizon/timestep);   % Number of discrete timesteps
time            = [0:timestep:(nSteps)*timestep];   % Create array with discrete time steps
nx              = 7;                                % Size of the state
nu              = 1;                                % Size of input u
foh             = 1;                                % First order hold reconstruction is on. (1)

weight = [1-Damagereduction Damagereduction];       % Weight for Weighted Sum MOOP
uMax            = 33^2;                             % input voltage is limited
% Construct the OCP
% The simulation is started with the flap already in motion. The initial
% state x0 represents the state after 100 seconds of swing in.
% x(1)   = first derivative of angle
% x(2)   = angle of the flap
% x(3-5) = is a state vector describing the waves generated by the device with its own motion
% x(6)   = Energy Harvested
% x(7)   = Damage aculmulated
x0 =   [0.0515     0.1392   314.4363   94.1062  190.4844         0         0 ]';

% import casadi and set solver
import casadi.*
ocp = casadi.Opti();
ocp.solver('ipopt')

% The wave differential equation
wave_dgl = @(x,u,d) [Ac * x(1:5) - Bc * 1e6 * u * gamma * x(2) + Bc * d;
    cost_energy(x,u,d);
    cost_damage(x,u);
    ];
% create the decision variables for the ocp.
x0_p = ocp.parameter(nx,1);
x = [x0_p, ocp.variable(nx, nSteps+1 - (foh))];
u = [ocp.parameter(nu, 1), ocp.variable(nu, nSteps)];
d = ocp.parameter(1, nSteps + foh);
h_integ = casadi.MX(nx, nSteps);

% integrate the discretised ode. x(t+1)  = integrator(x(t),u(t),d(t))
for iStep = 1:nSteps
    h_integ(:, iStep) = integrator_step_disturbed(x(:,iStep), u(:, iStep + (0:foh)), timestep, wave_dgl, d(:, iStep + (0:foh))) - x(:,iStep+1);
    ocp.subject_to( )
end
h_zero = h_integ(:) == 0;
ocp.subject_to( h_zero )

% Set initial values and apply wave distrurbance
try
    Wave_function = @(t) interp1(WaveData.time,WaveData.(['Wave_Seed_' num2str(Seed)]),t);
catch
    error('Please choose Seed between 1 and 10!')
end
ocp.set_value(d,arrayfun(@(t)Wave_function(t),time));
ocp.set_value(x0_p,x0)
ocp.set_value(u(1),0)
costfun = ([x(6,end), x(7,end)]);
p_params = ocp.parameter(2,1);
ocp.minimize( costfun*p_params );
ocp.set_value(p_params,weight);
ocp.subject_to(0<=u<=uMax);

% the MPC starts here. The MPC solutions are tracked in the
% ResultsMPC array. The applied signal is stored in the AppliedSignal array
% AppliedSignal (1,:) : Global time
% AppliedSignal (2,:) : Applied Voltage
% AppliedSignal (3,:) : Angle of Flap

starttime = [0:timestep:timehorizon];
ResultsMPC = cell(length(starttime),1);
AppliedSignal=[];
if (plotting)
    close all
    f = openfig("PlotBlank.fig");
    ax1 = f.Children(3);
    ax2 = f.Children(2);
    ax3 = f.Children(1);
end

for i = 1:length(starttime)
    nHorizon = round(MPCtimehorizon/timestep);     % Number of Pareto Points
    time = [starttime(i):timestep:(nHorizon)*timestep+starttime(i)];    % Create array with discrete time steps
    ocp.set_value(d,arrayfun(@(t) Wave_function(t),time));

    % if not the first MPC step set u_0 and x_0 to the solution of
    % previous MPC.
    if ~(i == 1)
        % reset solver for warmstarting

        ocp.set_value(x(:,1),solOld.value(x(:,2)));
        ocp.set_value(u(1),solOld.value(u(2)));

        ocp.set_initial(x(:,2:end-1),solOld.value(x(:,1:end-2)));
        ocp.set_initial(u(2:end-1),solOld.value(u(1:end-2)));

    end
    if i == 2
        % Set the solver to warmstart
        %  (faster because we have good initial guess from precious MPC)
        solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
            'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
            'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
        ocp.solver('ipopt', struct(), solver_opt)
    end
    try
        solOld = ocp.solve();
    catch
        % Faster solver sometimes cannot find solution.If that is the case
        % disable warmstart for one iteration
        ocp.solver('ipopt')
        solOld = ocp.solve();
        solver_opt = struct('warm_start_init_point', 'yes', 'mu_init', 1e-6, 'warm_start_mult_bound_push',1e-8, ...
            'warm_start_slack_bound_push', 1e-8, 'warm_start_bound_push', 1e-6, ...
            'warm_start_bound_frac',1e-6, 'warm_start_slack_bound_frac',1e-8);
        ocp.solver('ipopt', struct(), solver_opt)

    end

    %Save results in Array
    solMpc = struct;
    solMpc.x = solOld.value(x);
    solMpc.u = solOld.value(u);
    solMpc.time = solOld.value(time);
    AppliedSignal = [AppliedSignal [time(1); solOld.value(u(1)) ; solOld.value(x(2,1))]];
    ResultsMPC{i} = solMpc;


    if (plotting)



        axes(ax1)
        hold on

        plot(AppliedSignal(1,:),rad2deg(AppliedSignal(3,:)),'b',LineWidth=8)

        plot(solOld.value(time),rad2deg(solOld.value(x(2,:))),'r--',LineWidth=8)
        if ~(i==1)
            delete(ax1.Children(end))
            delete(ax1.Children(end))
        end

        axes(ax2)
        hold on

        plot(AppliedSignal(1,:),AppliedSignal(2,:),'b',LineWidth=8)

        plot(solOld.value(time),solOld.value(u),'r--',LineWidth=8)
        if ~(i==1)
            delete(ax2.Children(end))
            delete(ax2.Children(end))
        end

        axes(ax3)
        hold on

        plot(AppliedSignal(1,:),arrayfun(@(t) Wave_function(t),AppliedSignal(1,:)),'b',LineWidth=8)
        plot(solOld.value(time),arrayfun(@(t) Wave_function(t),solOld.value(time)),'r--',LineWidth=8)
        if ~(i==1)
            delete(ax3.Children(end))
            delete(ax3.Children(end))
        end
    end

end



function cost = cost_damage(x, u, ~, ~)
%     cost = (max(u - 484/(cos(x(2)).^2), 0).^2)*1e-6;
cost = (max(u - 484, 0).^2)*1e-6;
end

function cost = cost_energy(x, u, d)
persistent Ch S R0

if isempty(Ch) || isempty(S) || isempty(R0)
    load("PolySurge_inputs.mat" ,'Ch', 'S', 'R0');
end
cost = (Ch*x(1).^2 + x(3:5)'*S*x(3:5) - d .* x(1))*1e-6 + u/R0;
end



function x_end = integrator_step_disturbed(x0, u, dt, odefun, d)
% calculate one integration step with step size dt
import casadi.*

x0_rk = x0;
k = casadi.MX( size( x0, 1 ), 4 );

if size(u,2) == 1
    k(:,1) = odefun(x0_rk(:,end)                  , u, d);
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u, d);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u, d);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u, d);
else
    k(:,1) = odefun(x0_rk(:,end)                  , u(:,1),       d(:,1));
    k(:,2) = odefun(x0_rk(:,end) + dt / 2 * k(:,1), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,3) = odefun(x0_rk(:,end) + dt / 2 * k(:,2), u*[0.5; 0.5], d*[0.5; 0.5]);
    k(:,4) = odefun(x0_rk(:,end) + dt     * k(:,3), u(:,2),       d(:,2));
end

x_end  = x0_rk + dt / 6 * k * [1 2 2 1]';

end

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