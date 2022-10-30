function [Value] = FBM(time,args)
% This function creates and stores a simplex image
% Then it also handles the interpolation between
% 2 steps if "time" is not a descrete integer
arguments
    time                 (1,1) {mustBeNumeric}
    args.Seed            (1,1) {mustBeNumeric} = 1
    %Seed for the generation of the Simplex
    args.FBMSize         (1,1) {mustBeNumeric} = 999
    %Defines the size of the simplex input image
    args.warpFactor      (1,1) {mustBeNumeric} = 2
    %Scales the Distribution so the Wave can be more or less smooth
    args.predefinedMaps  (1,1) {boolean} = false;
end
if ~(args.predefinedMaps) 
persistent Map

    varmap = strcat('v',num2str(args.Seed));
    
    if ~isfield(Map, varmap)
        disp('New Map Generated')
        Map.(varmap) = FBMNoise(args.FBMSize,'Seed',args.Seed);
    end

FBMMap = Map.(varmap);
else
    % Currently 200 Maps Stored so seed is modded with 200
    try
    eval(['load PredefMap200.mat v' num2str(mod(args.Seed,200)) ';'])
    eval(['FBMMap = v' num2str(mod(args.Seed,200)) ';'])
    catch 
        warning('Run "RuntoCReatePredefinedMaps.m" first!')
    end 
       
end
time = time*args.warpFactor;  % Warp Factor is applied High warp means faster rng

%% Check if time is an integer and if not interpolate between the next two integers
if floor(time) == time
    Value = HamiltonianPathForFBMNew(FBMMap,time);
else
    ValueUpper = HamiltonianPathForFBMNew(FBMMap,ceil(time));
    ValueLower = HamiltonianPathForFBMNew(FBMMap,floor(time));
    Value = interp1([floor(time) ceil(time)],[ValueLower ValueUpper],time);
end
end

