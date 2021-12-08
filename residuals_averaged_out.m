function [system] = residuals_averaged_out(n_elements_k,n_elements_z, k_elements, ...
    Z, l_min, l_max, alpha, delta, k_min, k_max, beta, transition_matrix, ...
    theta, n_evaluations)
%residuals_averaged_out 
    %Inputs 
        % Parameters, number of elements 

 
system = zeros(n_elements_k * n_elements_z, 1); 

    for zz = 1:length(Z)

        equation = zeros(n_elements_k, 1);
        
        for kk = 1:(n_elements_k)
        
            points = transpose(linspace(k_elements(kk), k_elements(kk+1), n_evaluations));
            width = points(2) - points(1);
            fun_height_1 = zeros(n_evaluations,1); 
            fun_height_2 = zeros(n_evaluations,1);
        
            for ii = 1:(length(points))
                
                weight = weight_fun_capital(points(ii), k_elements);
                w_res = residual_fun(points(ii), Z(zz), theta, l_min, l_max, ...
                        alpha, delta, k_min, k_max, Z, ...
                        transition_matrix, k_elements, beta); 
                fun_height_1(ii,1) = weight(kk) * w_res;
                fun_height_2(ii,1) = weight(kk) * w_res;
        
            end
        
            equation(kk) = (sum(fun_height_1) + sum(fun_height_2)) * width;
        
        end
        
        if zz < 2
        
            system(1:8) = equation; 

        else
            
            if zz > 2

                system(17:24) = equation; 

            else

                system(9:16) = equation; 

            end

        end 
    
    end

end 



