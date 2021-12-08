function [res] = residual_fun(k, z, theta, l_min, l_max, alpha, ...
    delta, k_min, k_max, Z_vector, transition_matrix, k_elements, beta)
% weighted_residual Function
%   Give me states (k,z) and a vector theta_temp 
%   I give you the weighted residual 

% k has to make sense
if k > k_max 
    
    disp('k greater than the max value')
    stop

end

if k < k_min 

    dips('k smaller than the min value')
    stop

end

% Select the values of theta_temp that make sense 
theta_1 = theta(1:8,1);
theta_2 = theta(9:16,1);
theta_3 = theta(17:24,1); 
theta_temp = [theta_1, theta_2, theta_3]; 

% Position of state z
[~, z_pos] = min(abs(Z_vector - z)); 

% Labor today 
l = min(l_max,max(l_min, transpose(theta_temp(:,z_pos)) *  ...
                weight_fun_capital(k, k_elements)));

% Output today
y = exp(z) * k^(alpha) * l^(1-alpha);

% Consumption today
c = (1-alpha) * y / (l^2);

% Capital tomorrow 
kp = min(k_max, max(k_min, y + (1-delta) * k - c));  

% Recompute consumption today
c = max(10^(-7), y + (1-delta) * k - kp); 

% Recompute capital tomorrow 
kp = min(kp, y + (1-delta) * k); 

% This loop is over the three possible states tomorrow
lp = zeros(length(Z_vector),1);
yp = zeros(length(Z_vector),1);
cp = zeros(length(Z_vector),1);
kpp = zeros(length(Z_vector),1);

for zz = 1:(length(Z_vector))
    
    % Labor tomorrow 
    lp(zz,1) = min(l_max, max(l_min,transpose(theta_temp(:,zz)) *  ...
        weight_fun_capital(kp, k_elements)));

    % Output tomorrow
    yp(zz,1) = exp(zz) * kp^(alpha) * lp(zz,1)^(1-alpha);

    % Consumption tomorrow 
    cp(zz,1) = (1-alpha) * yp(zz,1) / (lp(zz,1))^2;

    % Capital in two periods 
    kpp(zz,1) = min(k_max, max(k_min, yp(zz,1) + (1-delta) * kp - cp(zz,1)));

    % Recompute consumption tomorrow 
    cp(zz,1) = max(10^(-7), yp(zz,1) + (1-delta) * kp - kpp(zz,1)); 

    % Recompute capital in two periods 
    kpp(zz,1) = min(kpp(zz,1), yp(zz,1) + (1-delta) * kp); 

end

% Now construct the residual 
res = beta * c * transition_matrix(z_pos, :)...
                * (1./cp .* (alpha * (yp / kp) + 1 - delta)) - 1;

% Choose the appropiate weight 
% [~, i] = min(abs(k-k_elements));
% k_elements_temp = k_elements; 
% k_elements_temp(i,1) = 1000000;
% [~, j] = min(abs(k-k_elements_temp));
% i = min(i,j);


% The value of k is between elements i and i+1
% weight = weight_fun_capital(k, k_elements);

% Vector of weighted residuals 
% w_res = weight * res; 

%w_res = weight(i) * res; 

%if i == 1

    %w_res_2 = 0; 

%else

    %w_res_2 = weight(i - 1) * res; 

%end



end