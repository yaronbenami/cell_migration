clear
clc

Da = 0.5;

phi_c_0 = [0.2 0.4]; % initial guess for the range in which phi_c is

options = optimset('TolX',1e-3);
phi_c = fzero(@(phi_c) calc_up_down_migration_prob_model(phi_c,Da), phi_c_0,options); 

fprintf('phi_c = %1.3f\n', phi_c);


