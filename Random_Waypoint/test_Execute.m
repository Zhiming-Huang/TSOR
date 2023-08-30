%Testing Random Waypoint mobility model.
clear all;clc;close all;

s_input = struct('V_POSITION_X_INTERVAL',[10 30],...%(m)
                 'V_POSITION_Y_INTERVAL',[10 30],...%(m)
                 'V_SPEED_INTERVAL',[3 3],...%(m/s)
                 'V_PAUSE_INTERVAL',[0 1],...%pause time (s)
                 'V_WALK_INTERVAL',[4.00 6.00],...%walk time (s)
                 'V_DIRECTION_INTERVAL',[-180 180],...%(degrees)
                 'SIMULATION_TIME',4000,...%(s)
                 'NB_NODES',4);
s_mobility = Generate_Mobility(s_input);

timeStep = 1;%(s)
% test_Animate(s_mobility,s_input,timeStep);
communication_Radius = [10,4,6,8,10,10];
% generate_Obs(s_mobility,s_input,timeStep,communication_Radius);
generate_Obs_fixed(s_mobility,s_input,timeStep,communication_Radius);