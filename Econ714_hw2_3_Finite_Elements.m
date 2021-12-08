% -------------------------------------------------------------------------
% Econ 714
% Homework 2
% Fall 2021 
% Javier Tasso
% Finite Elements
% -------------------------------------------------------------------------
clearvars
clc 
cd 'C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2'
% -------------------------------------------------------------------------

% Set values of parameters 
alpha = 0.33;
beta = 0.97;
delta = 0.1;

% Compute the steady state of the economy
    options = optimset('Display','off');
    F = @(x) [(1-alpha) / max(0, x(1)^(alpha) * x(2)^(1-alpha) - delta * x(1)) ...
        * (x(1)/x(2))^(alpha) - x(2); ... 
        (1/(alpha*beta) - (1-delta)/alpha) * (x(1)/x(2))^(1-alpha) - 1];
    [x,~] = fsolve(F,[1;1],options);

    % Check that everything looks ok
    if x(1)<0 || x(2)<0 || x(2)>1
        
        disp('There may be corner solutions, the program stops here')
        stop

    end

    % Store the solution 
    cons_ss = max(0, x(1)^(alpha) * x(2)^(1-alpha) - delta * x(1)); 
    kapi_ss = x(1);
    labo_ss = x(2); 
    prod_ss = kapi_ss^(alpha) * labo_ss^(1-alpha); 
    valu_ss = (log(cons_ss) - (labo_ss)^2 / 2);

    % Check 
    if abs(prod_ss - (cons_ss + delta * kapi_ss)) > 10^(-7)

        disp('There is something wrong with the steady state values, the program stops here')
        stop 

    end

% Clean and keep the steady state 
clear F prod_ss options x 

% -------------------------------------------------------------------------
% Define boundaries for capital 
k_min = 0.7 * kapi_ss;
k_max = kapi_ss * 1.3;
amplitude = k_max - k_min; 

% Define boundaries for labor 
l_min = 0.5;
l_max = 2;

% Productivity vector and its transition matrix 
z = [-0.05, 0, 0.05];
transition_matrix = [0.97, 0.03, 0; 0.01, 0.98, 0.01; 0, 0.03, 0.97];

% Define the finite elements for capital 
    % Here I use the same values for the three leves of productivity 
k_0 = k_min;
k_1 = (kapi_ss - k_min) / 2 + k_0;
k_2 = (kapi_ss - k_min) / 4 + k_1; 
k_3 = (kapi_ss - k_min) / 8 + k_2; 
k_4 = kapi_ss;
k_8 = k_max;
k_7 = k_8 - (k_max - kapi_ss) / 2; 
k_6 = k_7 - (k_max - kapi_ss) / 4;
k_5 = k_6 - (k_max - kapi_ss) / 8; 
k_elements = transpose([k_0, k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8]); 
clear k_0 k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 

% Need to define 8 psy functions
% See function weight_fun_capital(k, k_elements)

% Number of elements per productivity level 
n_elements_k = length(k_elements) - 1;
n_elements_z = length(z);
n_elements = n_elements_k * n_elements_z;

% Check that the weigth function makes sense
% k_grid = transpose(linspace(k_min,k_max, 1000));
% phi = zeros(1000, 8);
% for ii = 1:1000

    % phi(ii,:) = transpose(weight_fun_capital(k_grid(ii), k_elements)); 

 % end

%plot(k_grid, phi(:,1))
%hold on

%for tt = 2:8

    %plot(k_grid, phi(:,tt))

%end
%hold off 




[res] = residual_fun(4.5,-0.05, rand(24,1),0.5, 1.3, alpha, delta, k_min,k_max, z, ...
     transition_matrix, k_elements, beta);

[res2] = residual_fun(4.5,0, rand(24,1),0.5, 1.3, alpha, delta, k_min,k_max, z, ...
     transition_matrix, k_elements, beta);

[res3] = residual_fun(4.5,0.05, rand(24,1),0.5, 1.3, alpha, delta, k_min,k_max, z, ...
     transition_matrix, k_elements, beta);



 

% disp(wres)




 
% For a given theta evaluate the residual 
% [aaa] = residuals_averaged_out(n_elements_k, n_elements_z, k_elements, z, l_min, l_max, ...
    % alpha, delta, k_min, k_max, beta, transition_matrix, zeros(24), 100); 
%disp(aaa)

% theta_prev = zeros(24,1); 
% val_obj_fun = norm(residuals_averaged_out(n_elements_k, n_elements_z, k_elements, z, l_min, l_max, ...
    % alpha, delta, k_min, k_max, beta, transition_matrix, theta_prev, 10)); 
% val_obj_fun_prev = val_obj_fun; 
% n_it = 0;
 
 
% while -val_obj_fun < -0.01 && n_it < 10000
    
    % Draw new value of theta
    % theta_temp = theta_prev + normrnd(0,1,24,1);

    % Compute value at this new value of theta
    % val_temp = norm(residuals_averaged_out(n_elements_k, n_elements_z, k_elements, z, l_min, l_max, ...
    % alpha, delta, k_min, k_max, beta, transition_matrix, theta_temp, 10)); 

    % Probability
    % p = min(1, val_temp / val_obj_fun_prev); 

    % Update the value of theta
    % bino_rand = binornd(1, p, 1);  
    % theta = theta_temp * bino_rand + theta_prev * (1 - bino_rand); 

     

    % update the value of the function if necessary 
    % if theta == theta_temp 

        % val_obj_fun = val_temp; 

    %end 
    


    % Update things 
    % n_it = n_it + 1; 
    % val_obj_fun_prev = val_temp; 
    % theta_prev = theta; 





% end

% stop 




 
tol_level = 10^(-3);
options=optimset('Display','Iter','TolFun',tol_level,'TolX',tol_level,...
    'MaxFunEvals', 10^5, 'MaxIter', 10^4);
theta = fsolve(@(x) residuals_averaged_out(n_elements_k, n_elements_z, ...
    k_elements, z, l_min, l_max, alpha, delta, k_min, k_max, beta, ...
    transition_matrix, x, 10), ones(24,1), options); 


% Maybe try to solve it using something random 

%[bbb] = residuals_averaged_out(n_elements_k, n_elements_z, k_elements, z, l_min, l_max, ...
    %alpha, delta, k_min, k_max, beta, transition_matrix, theta, 10); 

%disp(aaa)
 
% Suppose it works and continue 
k_grid = transpose(linspace(k_min,k_max, 1000)); 
L = zeros(length(k_grid),3); 

for kk = 1:length(k_grid)

    L(kk, 1) = transpose(theta(1:8)) * weight_fun_capital(k_grid(kk), k_elements); 
    L(kk, 2) = transpose(theta(9:16)) * weight_fun_capital(k_grid(kk), k_elements); 
    L(kk, 3) = transpose(theta(17:24)) * weight_fun_capital(k_grid(kk), k_elements); 

end

Y = k_grid.^(alpha) .* L.^(1-alpha); 
C = (1-alpha) * Y ./ L; 
Kp = min(k_max, max(k_min, Y + (1-delta) * k_grid - C)); 
C = max(10^(-7), Y + (1-delta) * Kp - Kp); 
V = (log(C) - (L.^2)./2); 

k_grid_plot = k_grid(400:600);
V_plot = V(400:600,:);

% Value function 
figure(1)
plot(k_grid_plot, V_plot(:,1), 'linewidth', 2)
hold on
plot(k_grid_plot, V_plot(:,2), 'linewidth', 2)
plot(k_grid_plot, V_plot(:,3), 'linewidth', 2)
xlabel('Capital')
ylabel('Value')
title('Value Functions')
hold off
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question6_plot_value_functions.png');
close(figure(1))


clearvars 




%%
stop 

% This does not work. I don't know what is happening 

 

integrals = zeros(n_elements_k, 3); 

for zz = 1:n_elements_z 

    for jj = 1:n_elements_k
    
        points = transpose(linspace(k_elements(jj), k_elements(jj+1), 1000));
        width = points(2) - points(1);
        fun_height = zeros(1000,1); 
        
        for ii = 1:length(points)
        
            fun_height(ii) = residual_fun(points(ii), z(zz), zeros(24,1), l_min, l_max, ...
                alpha, delta, k_min, k_max, z, ...
                transition_matrix, k_elements, beta); 
        
        end
        
        integrals(jj,zz) = sum(fun_height) * width; 
    
    end

end

clear fun_height jj ii points width 
integrals = [integrals(:,1); integrals(:,2); integrals(:,3)]; 
disp(integrals)
 
  









