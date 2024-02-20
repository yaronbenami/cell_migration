clc
clear

%% parameters
phi_ref = 0.4;
K = 10;
Pe_c = 300;
M = 1*K;
Da = 0.5;
T = 10;
U_c = 0.003;
Pe = 1;


%% x,t, and differentiation matrixes
x_max = 12;

N_x = 1000;
x = linspace(-x_max, x_max, N_x)';
dx = x(2) - x(1);

dt = 1*dx/U_c;
t = 0:dt:1000;

D_1 = (-0.5*diag(ones(1,N_x-1), -1) + 0.5*diag(ones(1,N_x-1), 1))/dx;
D_1(1,1:3) = [-1.5 2 -0.5]/dx;
D_1(end,end-2:end) = [0.5 -2 1.5]/dx;

D_2 = (1*diag(ones(1,N_x-1), -1) - 2*diag(ones(1,N_x), 0) + 1*diag(ones(1,N_x-1), 1))/dx^2;
D_2(1,1:4) = [2 -5 4 -1]/dx^2;
D_2(end,end-3:end) = [-1 4 -5 2]/dx^2;

O = zeros(N_x);
I = eye(N_x);


%% variables
phi = zeros(length(x),length(t));  % cell volume fraction
c = zeros(length(x),length(t));    % chemoattractnat concentration
psi = zeros(length(x),length(t));  % cell flux
phi_dif = zeros(length(x),length(t));  % difference between downstream and upstream volume fraction
psi_dif = zeros(length(x),length(t));  % difference between downstream and upstream flux

%% initial state
phi(:,1) = phi_ref*exp(-x.^2);
fun = @(z,Da,phi) exp(-z.^2+Da*phi*erf(z)*sqrt(pi)/2);
for jj = 1:length(x)
    c(jj,1) = Da*phi_ref*exp(-Da*phi_ref*sqrt(pi)/2*erf(x(jj)))*integral(@(z) fun(z,Da,phi_ref),-Inf,x(jj));
end


%% time loop

for ii = 1:length(t)-1

    phi_p = phi(:,ii);
    c_p = c(:,ii);
    psi_p = psi(:,ii);
    phi_dif_p = phi_dif(:,ii);
    psi_dif_p = psi_dif(:,ii);

    s_p = -K*phi_p.^2./(1-phi_p).^4 + M*D_1*c_p;
    F_p = coth(s_p)-1./s_p;
    F_p(abs(s_p) < 1e-4) = 1/3*s_p(abs(s_p) < 1e-4);  % limit - prevent singularity
    S_dif_p = (1 - exp(-s_p))./(1 + exp(-s_p));
    F_dif_p = (exp(-s_p).*s_p + exp(-s_p) + s_p - 1)./(s_p.*(exp(-s_p) + 1));
    F_dif_p(abs(s_p) < 1e-4) = 1/2 + 1/24*(s_p(abs(s_p) < 1e-4)).^2;   % limit - prevent singularity

    Eq_Mat = [I/dt - 1/Pe_c*D_2, O, D_1 
              -Da*I, I/dt + D_1 + Da*diag(phi_p), O
              -U_c/T*diag(F_p) , O, I*(1/dt + 1/T) - 1/Pe_c*D_2];

    RHS = [phi_p/dt  
           c_p/dt
           psi_p/dt];

    % BCs
    % no diffusive flux of cells
    Eq_Mat(1,:) = [D_1(1,:), O(1,:), O(1,:)];
    RHS(1) = 0;
    Eq_Mat(N_x,:) = [D_1(N_x,:), O(1,:), O(1,:)];
    RHS(N_x) = 0;
    % no chemokine inlet/outlet flux
    Eq_Mat(N_x+1,:) = [O(1,:), D_1(1,:), O(1,:)];
    RHS(N_x+1) = 0;
    % no advective flux of cells
    Eq_Mat(2*N_x+1,:) = [O(1,:), O(1,:), I(1,:)];
    RHS(2*N_x+1) = 0;
    Eq_Mat(3*N_x,:) = [O(N_x,:), O(N_x,:), I(N_x,:)];
    RHS(3*N_x) = 0;


    W = Eq_Mat\RHS;


    phi(:,ii+1) = W(1:N_x);
    c(:,ii+1) = W(N_x+1:2*N_x);
    psi(:,ii+1) = W(2*N_x+1:3*N_x);

    Eq_Mat_dif = [I*(1/dt + 1/T) - 1/Pe_c*D_2, D_1  
                   O , I*(1/dt + 1/T) - 1/Pe_c*D_2];

    RHS_dif = [phi_dif_p/dt + 1/T*phi_p.*S_dif_p 
               psi_dif_p/dt + U_c/T*phi_p.*F_dif_p];

    % BCs
    % no diffusive flux of cells
    Eq_Mat_dif(1,:) = [I(1,:), O(1,:)];
    RHS_dif(1) = 0;
    Eq_Mat_dif(N_x,:) = [I(N_x,:), O(1,:)];
    RHS(N_x) = 0;
    % no advective flux of cells
    Eq_Mat_dif(N_x+1,:) = [O(1,:), I(1,:)];
    RHS_dif(N_x+1) = 0;
    Eq_Mat_dif(2*N_x,:) = [O(N_x,:), I(N_x,:)];
    RHS_dif(2*N_x) = 0;

    W_dif = Eq_Mat_dif\RHS_dif;
    phi_dif(:,ii+1) = W_dif(1:N_x);
    psi_dif(:,ii+1) = W_dif(N_x+1:2*N_x);


    ii

end

% plot at t = 0, 100, 500, 1000

[gar, ind_1] = min(abs(t - 100));
[gar, ind_2] = min(abs(t - 500));
ind_t = [1,ind_1, ind_2, length(t)]; %round(length(t)/6*(1:6));


figure(1)
plot(x, phi(:,ind_t))
xlabel('x')
ylabel('\phi')

figure(2)
plot(x, psi(:,ind_t)./phi(:,ind_t)/U_c)
xlabel('x')
ylabel('u_c')

figure(3)
plot(x, c(:,ind_t))
xlabel('x')
ylabel('c')

















