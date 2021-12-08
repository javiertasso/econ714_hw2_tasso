% -------------------------------------------------------------------------
% Econ 714
% Homework 2
% Fall 2021 
% Javier Tasso
% Chebyshev Polynomials
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

% Define boundaries for capital 
cover_grid = 0.3;
k_min = kapi_ss * (1 - cover_grid);
k_max = kapi_ss * (1 + cover_grid);
interval = kapi_ss * 2 * cover_grid; 
l_min = labo_ss * (1-cover_grid);
l_max = labo_ss * (1 + cover_grid);

% Write down the transition matrix
z = [-0.05, 0, 0.05];
transition_matrix = [0.97, 0.03, 0; 0.01, 0.98, 0.01; 0, 0.03, 0.97];

% Set up the number of polynomials 
n_polynomials = 5;
shock_num = length(z);  

% Set timer
tic

% Get the order of polynomials
M = n_polynomials * shock_num; 

% Find zeros of the Chebychef Polynomial of order M
ZC = -cos((2*(1:n_polynomials)'-1) * pi / (2 * n_polynomials)); 

% Define Chebychef Polinomials
T_k=ones(n_polynomials,n_polynomials);
T_k(:,2)=ZC;

for ii=3:n_polynomials

        T_k(:,ii)=2*ZC.*T_k(:,ii-1)-T_k(:,ii-2);

end

% Project collocation points in the K space
grid_k=((ZC+1)*(k_max-k_min))/2+k_min;

% Initial Guess for Chebyshev coefficients
rho_guess = zeros(2*M,1);

for zz = 1:shock_num

    rho_guess((zz-1) * n_polynomials + 1) = valu_ss;
    rho_guess((zz-1) * n_polynomials + 1 + M) = labo_ss; 
   
end

 
  

% Solve for Chebyshef coefficients
% Here we use the residual_fcn
tol_level = 10^(-10);
rho = residual_fcn_2(alpha,beta,delta,k_min,k_max,rho_guess,...
    grid_k,T_k,z,transition_matrix,n_polynomials,shock_num,M, l_min, l_max, tol_level);
rho_old = rho;
n_polynomials_old = n_polynomials;
time_cheby = toc; 



% Once we have the coefficients compute the value functions. I have to do
% this later because we need function eulerr_grid, function simulation 
grid_num = 1000; 
grid_k_complete = zeros(grid_num,1);
for i = 1:grid_num
    grid_k_complete(i) = k_min + (i-1)*interval/(grid_num-1);
end

  

[g_k,g_c,g_l,value_fcn,euler_error,max_error]= ...
                    eulerr_grid_2(alpha,beta,delta,rho,z,transition_matrix,...
                    k_min,k_max,grid_k_complete,shock_num,n_polynomials,grid_num,M);

figure(1)
plot(grid_k_complete,euler_error)
title('Log10 Euler Error')

% Euler Error
figure(1)
plot(grid_k_complete,euler_error(:,1), 'linewidth', 2)
hold on
plot(grid_k_complete,euler_error(:,2), 'linewidth', 2)
plot(grid_k_complete,euler_error(:,3), 'linewidth', 2)
title('Log10 Euler Error')
subtitle(['Time (seconds): ', num2str(time_cheby), ...
      ' - Tolerance level: ', num2str(tol_level)])
xlabel('Capital')
ylabel('Euler Equation Error')
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question5_plot_eee.png');
close(figure(1))

% Value function
figure(2)
plot(grid_k_complete,value_fcn(:,1), 'linewidth', 2)
hold on
plot(grid_k_complete,value_fcn(:,2), 'linewidth', 2)
plot(grid_k_complete,value_fcn(:,3), 'linewidth', 2)
title('Value Functions')
subtitle(['Time (seconds): ', num2str(time_cheby), ...
      ' - Tolerance level: ', num2str(tol_level)])
xlabel('Capital')
ylabel('Value')
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question5_plot_value_fun.png');
close(figure(2))

% Policy function
figure(3)
plot(grid_k_complete,g_k(:,1), 'linewidth', 2)
hold on
plot(grid_k_complete,g_k(:,2), 'linewidth', 2)
plot(grid_k_complete,g_k(:,3), 'linewidth', 2)
plot(grid_k_complete, grid_k_complete, 'linewidth', 1, 'color', 'black')
title('Policy Functions')
subtitle(['Time (seconds): ', num2str(time_cheby), ...
      ' - Tolerance level: ', num2str(tol_level)])
xlabel('Capital')
ylabel('Capital next period')
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question5_plot_pol_fun.png');
close(figure(3))



% Comparison 
clearvars -except value_fcn grid_k_complete
 

% Import data
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 13);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Kgrid", "V1", "V2", "V3", "K1", "K2", "K3", "L1", "L2", "L3", "n_iter", "time", "tol_level"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Import the data
    output_q2 = readtable("C:\Users\Javier\Dropbox\Aplicaciones\Overleaf\Econ714-PS2\econ714_homework2_question2_output.csv", opts);
    
    
    % Clear temporary variables
    clear opts

% Here make plots to compare with previos results 
% Compare with value function of question 2
% Focus on z=0

K_grid = output_q2 {:, 'Kgrid'};
% policy_fun_q2 = output_q2 {:, 'K2'};
% policy_fun_q2 = [K_grid, policy_fun_q2];
value_fun_q2 = output_q2 {:, 'V2'};
aaa = output_q2 {:, 'V1'}; 
value_fun_q2 = [aaa, value_fun_q2];
aaa = output_q2 {:, 'V3'};
value_fun_q2 = [value_fun_q2, aaa];
clear aaa
 
% value_fun_q2 = [K_grid, value_fun_q2];
% policy_fun_q5 = 

 
int_value_fun_q2 = zeros(length(grid_k_complete),3); 

for zz = 1:3

    for kk = 1:length(grid_k_complete)
    
        [~,i] = min(abs(grid_k_complete(kk) - K_grid));
        K_grid_temp = K_grid; 
        K_grid_temp(i) = 1000000;
        [~,j] = min(abs(grid_k_complete(kk) - K_grid_temp));
        
        i = min(i,j);
        int_value_fun_q2(kk,zz) = value_fun_q2(i,zz) + (value_fun_q2(i+1,zz) - value_fun_q2(i,zz)) / ...
            (K_grid(i+1) - K_grid(i)) * (grid_k_complete(kk) - K_grid(i));
             
    end

end


% Comparison with question 2
figure(4)
plot(grid_k_complete,abs(int_value_fun_q2(:,1) - value_fcn(:,1)), 'linewidth', 2)
hold on
plot(grid_k_complete,abs(int_value_fun_q2(:,2) - value_fcn(:,2)), 'linewidth', 2)
plot(grid_k_complete,abs(int_value_fun_q2(:,3) - value_fcn(:,3)), 'linewidth', 2)
title('Comparison with simple VFI')
% subtitle(['Time (seconds): ', num2str(time_cheby), ...
      % ' - Tolerance level: ', num2str(tol_level)])
xlabel('Capital')
ylabel('Difference (Absolute Value)')
legend({'z=-0.05','z=0','z=0.05'},'Location', 'southeast')
saveas(gcf,'econ714_homework2_question5_plot_comparison_q2.png');
close(figure(4))



clear i j kk K_grid_temp 


 

% Clean 
clearvars 

