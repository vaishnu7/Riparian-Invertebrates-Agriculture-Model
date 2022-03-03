clc;
clear;
close all;
format long
%%parameter values that will not vary
s=1; K=1; nu=1; L=1;
%%PARAMETERS WITH hopf at beta=67.988 approx.
% r=.002; alpha=.005; beta=67.9884746653243112746167753357440233230590820312503;
% theta=.9; gamma=.0005; delta=.003;
%%PARAMETERS WITH PERIODIC OSCILLATIONS
% r=.09; alpha=.005; beta=405.65645000000003504;
% theta=.9; gamma=.0005; delta=.003;
%%PARAMETERS WITH PERIODIC OSCILLATIONS
r=.9; alpha=.005; beta=1.5;
theta=.9; gamma=.0005; delta=.003;
% % %PARAMETERS WITH ALL CONDITIONS OF GLOBAL STABILITY
% r=10; alpha=0.9;
% theta=.5; gamma=3; delta=.1;
% beta=4;
%main model solving in RK4 method
riparian_dim=@(t,x)[r*x(1)*(1 - (x(1)/K)) - alpha*(x(1)^2)*x(3) - beta*x(2)*x(1); theta*beta*x(2)*x(1) - gamma*x(2)^2 - delta*x(2)*x(3); s*x(3)*(1 - (x(3)/L)) + nu*x(2)*x(3)];
[t,x]=ode45(riparian_dim,[0 1000],[0 0 0]);
p_riparian=size(x);
V_cap_rev=x(p_riparian(1),1)
I_cap_rev=x(p_riparian(1),2)
A_cap_rev=x(p_riparian(1),3)
% %Equilibrium Existence
a_cap = (alpha*L*nu*theta*beta)/(s*(gamma+(delta*L*nu)/s))
b_cap = (r/K) + alpha*L -(alpha*(L^2)*nu*delta/(s*(gamma + (delta*L*nu/s)))) + (theta*(beta^2)/(gamma + (delta*L*nu/s)))
c_cap = -r - (delta*beta*L/(gamma + (delta*L*nu/s)))
D_cap = b_cap^2-4*a_cap*c_cap;
V1_cap = (-b_cap+sqrt(D_cap))/(2*a_cap)
V2_cap = (-b_cap-sqrt(D_cap))/(2*a_cap)
V_cap_act = V1_cap
I_cap_act = (theta*beta*V_cap_act - delta*L) / (gamma + (delta*L*nu/s))
A_cap_act = (s + nu*I_cap_act)*(L/s)
% %SENSITIVITY WITH RESPECT TO BETA
beta_sen=@(t,x_beta)[r*x_beta(1)-2*r*V_cap_act*x_beta(1)-alpha*(2*V_cap_act*A_cap_act*x_beta(1)+V_cap_act^2*x_beta(3))-I_cap_act*V_cap_act-beta*(x_beta(2)*V_cap_act+I_cap_act*x_beta(1));theta*I_cap_act*V_cap_act+theta*beta*(x_beta(2)*V_cap_act+I_cap_act*x_beta(1))-2*gamma*I_cap_act*x_beta(2)-delta*(x_beta(2)*A_cap_act+I_cap_act*x_beta(3));x_beta(3)-2*A_cap_act*x_beta(3)+x_beta(2)*A_cap_act+I_cap_act*x_beta(3)];
[t,x_beta]=ode45(beta_sen,[0 9000],[0 0 0])
x_beta_index=beta*x_beta
x_beta_log_index(:,1)=x_beta_index(:,1)/V_cap_act
x_beta_log_index(:,2)=x_beta_index(:,2)/I_cap_act
x_beta_log_index(:,3)=x_beta_index(:,3)/A_cap_act
figure(10)
plot(t,x_beta_index(:,1))
hold on
figure(11)
plot(t,x_beta_index(:,2))
hold on
figure(12)
plot(t,x_beta_index(:,3))
hold on
figure(14)
plot(t,x_beta_log_index(:,1))
hold on
figure(15)
plot(t,x_beta_log_index(:,2))
hold on
figure(16)
plot(t,x_beta_log_index(:,3))
hold on
q_beta=size(x_beta)
V_beta=x_beta(q_beta(1),1),I_beta=x_beta(q_beta(1),2),A_beta=x_beta(q_beta(1),3)
sen_index(1,1)=V_beta*alpha/V_cap_act;
sen_index(1,2)=I_beta*alpha/I_cap_act;
sen_index(1,3)=A_beta*alpha/A_cap_act;
sen_index(1,:)