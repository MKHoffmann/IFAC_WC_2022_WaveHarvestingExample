function [waveDisturbance varargout]= FBMStochasticWave(t, args)


arguments    
    t               (1,1) {mustBeNumeric}
    args.T_wave     (1,1) {mustBeNumeric} = 10       % Period of dominant Wave 
    args.H_wave     (1,1) {mustBeNumeric} = 0.35      % Hight factor
    args.TargetAmp  (1,1) {mustBeNumeric} = 4.5171e+06; %  From harmonic wave 
    args.N_freq     (1,1) {mustBeNumeric} = 50       %Number of Frequencies in Wave
    args.d_W        (1,1) {mustBeNumeric} = pi/20
    args.FreqRange  (2,1) {mustBeNumeric} = [0.05 3] %Range of frequencies in terms of multiple of T_wave
    args.AmpWarp    (1,1) {mustBeNumeric} = 2
    args.PhaseWarp  (1,1) {mustBeNumeric} = 0.2
    args.Seed       (1,1) {mustBeNumeric} = 1
    args.DrawDebug  (1,1) {boolean} = false;         %If selected the wave is drawn forever
    args.predefinedMaps  (1,1) {boolean} = false;    
end
persistent Gamma_nu
persistent Power
d_W = args.d_W;
N_freq = args.N_freq;  
H_wave = args.H_wave;
T_wave = args.T_wave;




w_max= args.FreqRange(2)*2*pi*(1/T_wave);  %% What Range of frequencies will be present in Wave
w_min= args.FreqRange(1);


% if isempty('Gamma_nu')
    load('PolySurge_inputs.mat', 'Gamma_nu'); 
    load('PolySurge_inputs.mat', 'Gamma_lu');
% end
%% Function for the amplitude of the wave depending on the frequency 
waveSpectrum = @(w) 262.9*H_wave^2*T_wave^(-4)*w^(-5)*exp(-1054*T_wave^(-4)*w^(-4));
%% Actual sample of wave frequencies
DataSample = (w_min:(w_max-w_min)/(N_freq-1):w_max);
%% Calculate AK
Ak =@(w) (2*d_W*waveSpectrum(w))^0.5;
% Calculate exitation coeficient 
% Gamma_nu = [[0  Gamma_nu(1,2)] ;Gamma_nu];
% Gamma = @(w) interp1(Gamma_nu(:,1), Gamma_nu(:,2), w);

p = polyfit(Gamma_lu.omega,Gamma_lu.Gamma,1);
Gamma = @(w) polyval(p,w);
 %% Rescale the Wave depending on the Sum of Amplitudes
varmap = strcat('f',num2str(args.N_freq));

if (~isfield(Power, varmap))
    Power.(varmap) = 1;
    disp('New Frequency scaled')
    Temp = zeros(100,1);
    for i = 1:1000
        
        [Temp(i),~] = FBMStochasticWave(i,'N_freq',N_freq,'Seed',args.Seed,'predefinedMaps',args.predefinedMaps);
    end
    Power.(varmap) = args.TargetAmp/mean(abs(Temp));
end

 %
waveDisturbance=0;
IndividualWaves = zeros(length(DataSample),1);

%% The final Wave is the sum over all the individual waves
%% Two sources of randomness are Amplitude offse and Phase Offset
%% Both calculated using an individual map and Seed
%%\\TODO Make persitent
seedPhrase = args.Seed;

rng(seedPhrase,'philox');
randomShift = rand(1,length(DataSample))*2*pi;

for i=1:length(DataSample)
%        AmplitudeOffset = 0.75+ 0.5*FBM(t,"Seed",(args.Seed+i+100),'warpFactor',args.AmpWarp);
        AmplitudeOffset = 1;
        PhaseOffset =    randomShift(i)+2*pi*FBM(t,"Seed",(args.Seed+i),'warpFactor',args.PhaseWarp,'predefinedMaps',args.predefinedMaps)  ;   
%         PhaseOffset =   2*pi*FBM(t,"Seed",(args.Seed+i+200),'warpFactor',args.PhaseWarp)  ;   
    
%            IndividualWaves(i)= Power.(varmap)*AmplitudeOffset*Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+PhaseOffset);
          IndividualWaves(i)= 1*AmplitudeOffset*Ak(DataSample(i)).*Gamma(DataSample(i)).*sin(DataSample(i).*t+PhaseOffset);

          %           IndividualWaves(i)=PhaseOffset;   

end
waveDisturbance = sum(IndividualWaves);
varargout= cell(3);
varargout{1} = arrayfun(@(w) Ak(w),DataSample).*arrayfun(@(w) Gamma(w),DataSample);
varargout{2} = IndividualWaves;
% varargout{3} =Prediction;

%% If Debug feature is selected the function enters into an infinite loop here
if args.DrawDebug
    increment = 0.5;
    f = figure(1);

    x1= [0 1];
    y1= [NewStochasticWave(0) NewStochasticWave(increment)];
    p = plot(x1,y1);
    p.XDataSource = 'x1';
    p.YDataSource = 'y1';
    k = 2*increment;
    
while(args.DrawDebug)
    [wave Individual]= SimplexStochasticWave(k,'N_freq',args.N_freq);
    x1 = [x1 k];
    y1 = [y1 wave];
    k = k+increment;
    refreshdata(f,'caller')
    drawnow
end
end

end