clear
clc

Da = 0.1;

phi_c_0 = [0.11 0.15]; % initial guess for phi_c

options = optimset('TolX',1e-3);
phi_c = fzero(@(phi_c) calc_up_down_migration_prob_model(phi_c,Da), phi_c_0,options); 

fprintf('phi_c = %1.3f\n', phi_c);

plot(Da,phi_c,'r*')

