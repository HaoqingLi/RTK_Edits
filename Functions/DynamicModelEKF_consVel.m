function [ out ] = DynamicModelEKF_consVel( dt, state )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = length(state);
out = state;
out(1:3) = out(1:3) + out(4:6)*dt;
end

