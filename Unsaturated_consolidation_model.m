clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approximate Analytical Solution

%% Foundamental parameters
g = 9.807; gamma_w = 9.807e3; T_at = 293.16; omega_a = 0.03;
u_atm = 1e5; beta_w = 4.2e-10; h = 0.02; R = 8.314;
m_a1 = 2e-7; m_a2 = -1e-7; m_w1 = 5e-8; m_w2 = 2e-7;
m_s1 = m_a1 + m_w1; m_s2 = m_a2 + m_w2; n_r = 0.5; S_r = 0.8;
k_ah = 1e-11; k_wh = 1e-10; k_as = 0.1*k_ah; k_ws = 0.1*k_wh;
k_av = 1e2*k_ah; k_wv = 1e2*k_wh; q = 1e5; H = 8;
r_w = 0.2; r_s = 0.6; r_e = 1.8; s = r_s/r_w; n = r_e/r_w;
%% Calculate the initial excess pore pressures
R_a1 = @(x)m_a2/(m_a2-m_a1-n_r*(1-S_r+h*S_r)/(x+u_atm));
R_a2 = @(x)m_a1/(m_a2-m_a1-n_r*(1-S_r+h*S_r)/(x+u_atm));
R_s1 = @(x)(m_s2-m_s1-n_r*(1-S_r+h*S_r)/(x+u_atm))/(m_s2+n_r*S_r*beta_w);
R_s2 = m_s1/(m_s2+n_r*S_r*beta_w);
f_a  = @(x)x-q*(R_s2*R_a1(x)-R_a2(x))/(1-R_s1(x)*R_a1(x));
u_a0 = fsolve(f_a,q/2);
u_w0 = q*(R_s2-R_s1(u_a0)*R_a2(u_a0))/(1-R_s1(u_a0)*R_a1(u_a0));
u_0  = [u_a0;u_w0];
%% Calculation parameters
% The elements of matrix Kh
C_ah = R*T_at*k_ah/(omega_a*g*((u_a0+u_atm)*(m_a1-m_a2)+n_r*(1-S_r)));
C_wh = k_wh/(gamma_w*m_w2);
C_as = R*T_at*k_as/(omega_a*g*((u_a0+u_atm)*(m_a1-m_a2)+n_r*(1-S_r)));
C_ws = k_ws/(gamma_w*m_w2);
% The elements of matrix C
C_wa = m_a2*(u_a0+u_atm)/((u_a0+u_atm)*(m_a1-m_a2)+n_r*(1-S_r));
C_aw = m_w1/m_w2-1;
% The elements of matrix N(z)
mu_a = @(z)(n*n/(n*n-1))*(log(n/s)+(k_ah/k_as)*log(s)-3/4) + ...
    (s*s/(n*n-1))*(1-s*s/(4*n*n))+((k_ah/k_as)/(n*n-1))*((s^4-1)/(4*n*n)-s*s+1) + ...
    (z*(2*H-z)/(r_w*r_w))*(k_ah/k_av)*(1-1/(n*n));
mu_w = @(z)(n*n/(n*n-1))*(log(n/s)+(k_wh/k_ws)*log(s)-3/4) + ...
    (s*s/(n*n-1))*(1-s*s/(4*n*n))+((k_wh/k_ws)/(n*n-1))*((s^4-1)/(4*n*n)-s*s+1) + ...
    (z*(2*H-z)/(r_w*r_w))*(k_wh/k_wv)*(1-1/(n*n));
% The elements of matrix A(z)
a_11 = @(z)-(2*C_ah/(r_e*r_e*mu_a(z)))/(C_aw*C_wa-1);
a_12 = @(z)(2*C_wa*C_wh/(r_e*r_e*mu_w(z)))/(C_aw*C_wa-1);
a_21 = @(z)(2*C_aw*C_ah/(r_e*r_e*mu_a(z)))/(C_aw*C_wa-1);
a_22 = @(z)-(2*C_wh/(r_e*r_e*mu_w(z)))/(C_aw*C_wa-1);
delta = @(z)(a_11(z)-a_22(z))*(a_11(z)-a_22(z))+4*a_12(z)*a_21(z);
% The two eigenvalues of the matrix A(z)
lambda_1 = @(z)(a_11(z)+a_22(z)-sqrt(delta(z)))/2;
lambda_2 = @(z)(a_11(z)+a_22(z)+sqrt(delta(z)))/2;
%% Calculate the excess pore pressures
u_aas_a = @(z,t)exp(-lambda_1(z)*t)*((lambda_2(z)-a_11(z))*u_a0 - ...
    a_12(z)*u_w0)/sqrt(delta(z)) - exp(-lambda_2(z)*t)*((lambda_1(z) - ...
    a_11(z))*u_a0-a_12(z)*u_w0)/sqrt(delta(z));
u_aas_w = @(z,t)exp(-lambda_1(z)*t)*((lambda_2(z)-a_22(z))*u_w0 - ...
    a_21(z)*u_a0)/sqrt(delta(z)) - exp(-lambda_2(z)*t)*((lambda_1(z) - ...
    a_22(z))*u_w0-a_21(z)*u_a0)/sqrt(delta(z));
%% Calculate the average degree of consolidation
% The details of average degree of consolidation calculated using 5 point
% Gaussian-Legendre quadrature method can refer to the website
% https://en.wikipedia.org/wiki/Gaussian_quadrature
% The nodes P_x and weights W_w are given below:
P_x = [-(1/3)*sqrt(5+2*sqrt(10/7));-(1/3)*sqrt(5-2*sqrt(10/7));0; ...
    (1/3)*sqrt(5-2*sqrt(10/7));(1/3)*sqrt(5+2*sqrt(10/7))];
W_w = [(322-13*sqrt(70))/900;(322+13*sqrt(70))/900;128/225; ...
    (322+13*sqrt(70))/900;(322-13*sqrt(70))/900];
% Average degree of consolidation
U_ave = @(t)1-((W_w(1)*((m_s1-m_s2)*u_aas_a((H/2)*(P_x(1)+1),t)+m_s2* ...
    u_aas_w((H/2)*(P_x(1)+1),t)))/(2*((m_s1-m_s2)*u_a0+m_s2*u_w0)) + ...
    (W_w(2)*((m_s1-m_s2)*u_aas_a((H/2)*(P_x(2)+1),t)+m_s2* ...
    u_aas_w((H/2)*(P_x(2)+1),t)))/(2*((m_s1-m_s2)*u_a0+m_s2*u_w0)) + ...
    (W_w(3)*((m_s1-m_s2)*u_aas_a((H/2)*(P_x(3)+1),t)+m_s2* ...
    u_aas_w((H/2)*(P_x(3)+1),t)))/(2*((m_s1-m_s2)*u_a0+m_s2*u_w0)) + ...
    (W_w(4)*((m_s1-m_s2)*u_aas_a((H/2)*(P_x(4)+1),t)+m_s2* ...
    u_aas_w((H/2)*(P_x(4)+1),t)))/(2*((m_s1-m_s2)*u_a0+m_s2*u_w0)) + ...
    (W_w(5)*((m_s1-m_s2)*u_aas_a((H/2)*(P_x(5)+1),t)+m_s2* ...
    u_aas_w((H/2)*(P_x(5)+1),t)))/(2*((m_s1-m_s2)*u_a0+m_s2*u_w0)));
%% Output data
% Time factor Th
T_aas = [1e-6;1.333e-6;1.666e-6;2e-6;2.5e-6;3e-6;4e-6;5e-6;6e-6;7e-6;8e-6;9e-6; ...
    1e-5;1.333e-5;1.666e-5;2e-5;2.5e-5;3e-5;4e-5;5e-5;6e-5;7e-5;8e-5;9e-5; ...
    1e-4;1.333e-4;1.666e-4;2e-4;2.5e-4;3e-4;4e-4;5e-4;6e-4;7e-4;8e-4;9e-4; ...
    1e-3;1.333e-3;1.666e-3;2e-3;2.5e-3;3e-3;4e-3;5e-3;6e-3;7e-3;8e-3;9e-3; ...
    1e-2;1.333e-2;1.666e-2;2e-2;2.5e-2;3e-2;4e-2;5e-2;6e-2;7e-2;8e-2;9e-2; ...
    1e-1;1.333e-1;1.666e-1;2e-1;2.5e-1;3e-1;4e-1;5e-1;6e-1;7e-1;8e-1;9e-1; ...
    1;1.333;1.666;2;2.5;3;4;5;6;7;8;9; ...
    1e1;1.333e1;1.666e1;2e1;2.5e1;3e1;4e1;5e1;6e1;7e1;8e1;9e1;1e2];
% Time t corresponding to Th
t_aas = (gamma_w*m_w1*4*r_e*r_e/k_wh)*T_aas;
% Dimensionless excess pore-air and -water pressures
uD_aas_a = zeros(size(T_aas,1),1);
uD_aas_w = zeros(size(T_aas,1),1);
U_aas = zeros(size(T_aas,1),1);
for i = 1:size(T_aas,1)
    uD_aas_a(i) = u_aas_a(H/2,t_aas(i))/u_a0;
    uD_aas_w(i) = u_aas_w(H/2,t_aas(i))/u_w0;
    U_aas(i) = U_ave(t_aas(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Finite Difference Solution (without using Hansbo's treatment)

%% Foundamental parameters
global E J
E = diag([1 1]); % Matrix I in the paper (Eq. 37)
Dz = 0.2; % Space increment in the z-direction
J = H/Dz; % Space steps
Dt = @(l)200*l-100; % Time increments
L = 10000; % Time steps
t_fds = zeros(L,1);
T_fds = zeros(L,1);
for i = 1:L
    t_fds(i) = 100*i*i;
    T_fds(i) = t_fds(i)/(gamma_w*m_w1*4*r_e*r_e/k_wh);
end
m_s = [m_s1-m_s2 m_s2]; % Vector ms
mu_aa = (n*n/(n*n-1))*(log(n/s)+(k_ah/k_as)*log(s)-3/4) + ...
    (s*s/(n*n-1))*(1-s*s/(4*n*n))+((k_ah/k_as)/(n*n-1))*((s^4-1)/(4*n*n)-s*s+1);
mu_ww = (n*n/(n*n-1))*(log(n/s)+(k_wh/k_ws)*log(s)-3/4) + ...
    (s*s/(n*n-1))*(1-s*s/(4*n*n))+((k_wh/k_ws)/(n*n-1))*((s^4-1)/(4*n*n)-s*s+1);
N = diag([mu_aa mu_ww]); % Matrix N'
K_h = diag([C_ah C_wh]);
K_s = diag([C_as C_ws]);
C = [1 C_wa; C_aw 1];
B = diag([k_av/k_as k_wv/k_ws]);
P = @(l)(Dt(l)/(Dz*Dz*(n*n-1)))*(C^(-1))*K_s*B;
Q = @(l)(2*Dt(l)/(r_e*r_e))*((N*C)^(-1))*K_h;
V = @(l)(Q(l)^(-1))*P(l)+P(l);
W = @(l)E+2*V(l);
%% Define structure matrix
u_1 = cell(J+1,L); % store u_1,j^l
u_ave = cell(J+1,L); % store {\bar u}_j^l
u_ave_int = cell(J+1,1); % store d
u_fds_a = zeros(L,1);
u_fds_w = zeros(L,1);
U_fds = zeros(L,1); % store average degree of consolidation
%% Calculate dimensionless excess pore-air and -water pressures and average degree of consolidation
% the first time step
u_ave0 = cell(J+1,1);
u_ave0{1} = [0;0];
for j = 2:J+1
    u_ave0{j} = u_0;
end
u_int = chasing_method(V(1),W(1),u_ave0);
for j = 1:J+1
    u_1{j,1} = u_int{j};
end
u_ave{1,1} = (E-P(1)*(V(1)^(-1)))*u_0;
for j = 2:J
    u_ave{j,1} = u_ave0{j}+P(1)*u_1{j-1,1}-2*P(1)*u_1{j,1}+P(1)*u_1{j+1,1};
end
u_ave{J+1,1} = u_ave0{J+1}+2*P(1)*u_1{J,1}-2*P(1)*u_1{J+1,1};
for j = 1:J
    U_fds(1) = U_fds(1)+Dz*m_s*(u_ave{j,1}+u_ave{j+1,1})/(2*H*m_s*u_0);
end
U_fds(1) = 1-U_fds(1);
u_fds_a(1) = u_ave{0.5*J+1,1}(1);
u_fds_w(1) = u_ave{0.5*J+1,1}(2);
% Time steps from 2 to L
for k = 2:L
    k
    u_ave_int{1} = [0;0];
    for j = 2:J+1
        u_ave_int{j} = u_ave{j,k-1};
    end
    u_int = chasing_method(V(k),W(k),u_ave_int);
    for j = 1:J+1
        u_1{j,k} = u_int{j};
    end
    
    u_ave{1,k} = (E-P(k)*(V(k)^(-1)))*u_ave{1,k-1};
    for j = 2:J
        u_ave{j,k} = u_ave{j,k-1}+P(k)*u_1{j-1,k}-2*P(k)*u_1{j,k}+P(k)*u_1{j+1,k};
    end
    u_ave{J+1,k} = u_ave{J+1,k-1}+2*P(k)*u_1{J,k}-2*P(k)*u_1{J+1,k};
    for j = 1:J
    U_fds(k) = U_fds(k)+Dz*m_s*(u_ave{j,k}+u_ave{j+1,k})/(2*H*m_s*u_0);
    end
    U_fds(k) = 1-U_fds(k);
    u_fds_a(k) = u_ave{0.5*J+1,k}(1);
    u_fds_w(k) = u_ave{0.5*J+1,k}(2);
end
uD_fds_a(:) = u_fds_a(:)/u_a0;
uD_fds_w(:) = u_fds_w(:)/u_w0;
%% Figures of dimensionless excess pore-air pressure, pore-water pressure and average degree of consolidation
subplot(2,2,1);
semilogx(T_aas,uD_aas_a,'r',T_fds,uD_fds_a,'k');
xlim([1e-6,1e2])
set(gca,'ytick',0:0.2:1);
xlabel('T_h');
ylabel('u_a');
title 'Excess pore-air pressure';

subplot(2,2,2);
semilogx(T_aas,uD_aas_w,'r',T_fds,uD_fds_w,'k');
xlim([1e-6,1e2])
xlabel('T_h');
ylabel('u_w');
title 'Excess pore-water pressure';

subplot(2,2,3);
semilogx(T_aas,U_aas,'r',T_fds,U_fds,'k');
xlim([1e-6,1e2])
set(gca,'ytick',0:0.2:1);
xlabel('T_h');
ylabel('U');
title 'Average degree of consolidation';



