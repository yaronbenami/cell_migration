function DIFF  = calc_up_down_migration_prob_model(phi_ref, Da)


fprintf('phi = %1.3f\n', phi_ref);

%% parameters
K_M = 1;
K = 0.1;
Pe_c = 300;
M = K/K_M;
T = 10;
U_c = 0.003;


%% x,t, and differentiation matrixes
x_max = 12;

N_x = 1000;
x = linspace(-x_max, x_max, N_x)';
dx = x(2) - x(1);

dt = 0.5*dx/U_c;
t = 0:dt:100;

D_1 = (-0.5*diag(ones(1,N_x-1), -1) + 0.5*diag(ones(1,N_x-1), 1))/dx;
D_1(1,1:3) = [-1.5 2 -0.5]/dx;
D_1(end,end-2:end) = [0.5 -2 1.5]/dx;

D_2 = (1*diag(ones(1,N_x-1), -1) - 2*diag(ones(1,N_x), 0) + 1*diag(ones(1,N_x-1), 1))/dx^2;
D_2(1,1:4) = [2 -5 4 -1]/dx^2;
D_2(end,end-3:end) = [-1 4 -5 2]/dx^2;

%% integration matrix
int_matrix = zeros(N_x);
for jj=2:length(x)
    int_matrix(jj,1:jj) = dx*[0.5,ones(1,jj-2),0.5];
end

O = zeros(N_x);
I = eye(N_x);


%% variables
phi = zeros(length(x),length(t));
c = zeros(length(x),length(t));
psi = zeros(length(x),length(t));
N_diff = zeros(1,length(t));


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

    s_p = -K*phi_p.^2./(1-phi_p).^4 + M*D_1*c_p;
    F_p = (exp(-2*s_p).*s_p + exp(-2*s_p) + s_p - 1)./(s_p.*(1 - exp(-2*s_p)));
    F_p(s_p < 1e-4) = 1/3*s_p(s_p < 1e-4);  % limit - prevent singularity

    Eq_Mat = [I/dt - 1/Pe_c*D_2, D_1 
              -U_c/T*diag(F_p) , I*(1/dt + 1/T) - 1/Pe_c*D_2];

    RHS = [phi_p/dt  
           psi_p/dt];

    % BCs
    % no diffusive flux of cells
    Eq_Mat(1,:) = [D_1(1,:), O(1,:)];
    RHS(1) = 0;
    Eq_Mat(N_x,:) = [D_1(N_x,:), O(1,:)];
    RHS(N_x) = 0;
    % no advective flux of cells
    Eq_Mat(N_x+1,:) = [O(1,:), I(1,:)];
    RHS(N_x+1) = 0;
    Eq_Mat(2*N_x,:) = [O(N_x,:), I(N_x,:)];
    RHS(2*N_x) = 0;


    W = Eq_Mat\RHS;

    phi(:,ii+1) = W(1:N_x);
    psi(:,ii+1) = W(N_x+1:2*N_x);
    c(:,ii+1) = 1 - exp(-Da*int_matrix*phi(:,ii+1));

    s = -K*phi(:,ii+1).^2./(1-phi(:,ii+1)).^4 + M*Da*phi(:,ii+1).*(1 - c(:,ii+1));

    R = trapz(x, phi(:,ii+1).*(1 - exp(-s))./(1 + exp(-s)));
    N_diff(ii+1) = (N_diff(ii)/dt + R/T)/(1/dt + 1/T);


end

DIFF = N_diff(end);

fprintf('N_diff = %1.9f\n', DIFF);

end


















