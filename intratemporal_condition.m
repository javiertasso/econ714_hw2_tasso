function [outputArg1] = intratemporal_condition(z,k,k_next,alpha, delta, l)

% Intratemporal Condition
%   This function takes advantage of the intratemporal condition of the
%   social planner. We will use this function to compute l for given values
%   of the current state (k,z) and the next value of capital k_next 

    outputArg1 = (exp(z) * (1-alpha) * (k / l)^(alpha)) / max(10^(-10), ... 
        exp(z) * k^(alpha) * l^(1-alpha) - (k_next - (1-delta) * k)) - l;

end