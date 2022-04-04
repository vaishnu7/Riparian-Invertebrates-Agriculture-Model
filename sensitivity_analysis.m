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
% r=.9; alpha=.005; beta=1.5;
% theta=.9; gamma=.0005; delta=.003;
% % %PARAMETERS WITH ALL CONDITIONS OF GLOBAL STABILITY
r=10; alpha=0.9;
theta=.5; gamma=3; delta=.1;
beta=4;
%main model solving in RK4 method
riparian_dim=@(t,x)[r*x(1)*(1 - (x(1)/K)) - alpha*(x(1)^2)*x(3) - beta*x(2)*x(1); theta*beta*x(2)*x(1) - gamma*x(2)^2 - delta*x(2)*x(3); s*x(3)*(1 - (x(3)/L)) + nu*x(2)*x(3)];
[t,x]=ode45(riparian_dim,[0 20],[.2 .3 .5]);
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
%SENSITIVITY WITH RESPECT TO r
r_sen=@(t,x_r)[(V_cap_act+r*x_r(1))-(V_cap_act^2+2*r*V_cap_act*x_r(1))-alpha*(2*V_cap_act*A_cap_act*x_r(1)+V_cap_act^2*x_r(3))-beta*(V_cap_act*x_r(2)+I_cap_act*x_r(1));theta*beta*(V_cap_act*x_r(2)+I_cap_act*x_r(1))-2*gamma*I_cap_act*x_r(2)-delta*(A_cap_act*x_r(2)+I_cap_act*x_r(3));x_r(3)-2*A_cap_act*x_r(3)+(A_cap_act*x_r(2)+I_cap_act*x_r(3))];
[t,x_r]=ode45(r_sen,[0 20],[0 0 0])
x_r_index=r*x_r
x_r_log_index(:,1)=x_r_index(:,1)/V_cap_act
x_r_log_index(:,2)=x_r_index(:,2)/I_cap_act
x_r_log_index(:,3)=x_r_index(:,3)/A_cap_act
% figure(10)
% plot(t,x_r_index(:,1))
% hold on
% figure(11)
% plot(t,x_r_index(:,2))
% hold on
% figure(12)
% plot(t,x_r_index(:,3))
% hold on
% figure(14)
% plot(t,x_r_log_index(:,1))
% hold on
% figure(15)
% plot(t,x_r_log_index(:,2))
% hold on
% figure(16)
% plot(t,x_r_log_index(:,3))
% hold on
q_r=size(x_r)
V_r=x_r(q_r(1),1),I_r=x_r(q_r(1),2),A_r=x_r(q_r(1),3)
sen_index(1,1)=V_r*r/V_cap_act;
sen_index(1,2)=I_r*r/I_cap_act;
sen_index(1,3)=A_r*r/A_cap_act;
sen_index(1,:)
%SENSITIVITY WITH RESPECT TO alpha
alpha_sen=@(t,x_alpha)[(r*x_alpha(1))-(2*r*V_cap_act*x_alpha(1))-(V_cap_act^2*A_cap_act)-(alpha*(2*V_cap_act*A_cap_act*x_alpha(1)+V_cap_act^2*x_alpha(3)))-(beta*(x_alpha(2)*V_cap_act+I_cap_act*x_alpha(1)));(theta*beta*(x_alpha(2)*V_cap_act+I_cap_act*x_alpha(1)))-(2*gamma*I_cap_act*x_alpha(2))-(delta*(x_alpha(2)*A_cap_act+I_cap_act*x_alpha(3)));x_alpha(3)-(2*A_cap_act*x_alpha(3))+(x_alpha(2)*A_cap_act)+(I_cap_act*x_alpha(3))];
[t,x_alpha]=ode45(alpha_sen,[0 20],[0 0 0])
x_alpha_index=alpha*x_alpha
x_alpha_log_index(:,1)=x_alpha_index(:,1)/V_cap_act
x_alpha_log_index(:,2)=x_alpha_index(:,2)/I_cap_act
x_alpha_log_index(:,3)=x_alpha_index(:,3)/A_cap_act
% figure(10)
% plot(t,x_alpha_index(:,1))
% hold on
% figure(11)
% plot(t,x_alpha_index(:,2))
% hold on
% figure(12)
% plot(t,x_alpha_index(:,3))
% hold on
% figure(14)
% plot(t,x_alpha_log_index(:,1))
% hold on
% figure(15)
% plot(t,x_alpha_log_index(:,2))
% hold on
% figure(16)
% plot(t,x_alpha_log_index(:,3))
% hold on
q_alpha=size(x_alpha)
V_alpha=x_alpha(q_alpha(1),1),I_alpha=x_alpha(q_alpha(1),2),A_alpha=x_alpha(q_alpha(1),3)
sen_index(2,1)=V_alpha*alpha/V_cap_act;
sen_index(2,2)=I_alpha*alpha/I_cap_act;
sen_index(2,3)=A_alpha*alpha/A_cap_act;
sen_index(2,:)
%SENSITIVITY WITH RESPECT TO BETA
beta_sen=@(t,x_beta)[r*x_beta(1)-2*r*V_cap_act*x_beta(1)-alpha*(2*V_cap_act*A_cap_act*x_beta(1)+V_cap_act^2*x_beta(3))-I_cap_act*V_cap_act-beta*(x_beta(2)*V_cap_act+I_cap_act*x_beta(1));theta*I_cap_act*V_cap_act+theta*beta*(x_beta(2)*V_cap_act+I_cap_act*x_beta(1))-2*gamma*I_cap_act*x_beta(2)-delta*(x_beta(2)*A_cap_act+I_cap_act*x_beta(3));x_beta(3)-2*A_cap_act*x_beta(3)+x_beta(2)*A_cap_act+I_cap_act*x_beta(3)];
[t,x_beta]=ode45(beta_sen,[0 20],[0 0 0])
x_beta_index=beta*x_beta
x_beta_log_index(:,1)=x_beta_index(:,1)/V_cap_act
x_beta_log_index(:,2)=x_beta_index(:,2)/I_cap_act
x_beta_log_index(:,3)=x_beta_index(:,3)/A_cap_act
% figure(10)
% plot(t,x_beta_index(:,1))
% hold on
% figure(11)
% plot(t,x_beta_index(:,2))
% hold on
% figure(12)
% plot(t,x_beta_index(:,3))
% hold on
% figure(14)
% plot(t,x_beta_log_index(:,1))
% hold on
% figure(15)
% plot(t,x_beta_log_index(:,2))
% hold on
% figure(16)
% plot(t,x_beta_log_index(:,3))
% hold on
q_beta=size(x_beta)
V_beta=x_beta(q_beta(1),1),I_beta=x_beta(q_beta(1),2),A_beta=x_beta(q_beta(1),3)
sen_index(3,1)=V_beta*alpha/V_cap_act;
sen_index(3,2)=I_beta*alpha/I_cap_act;
sen_index(3,3)=A_beta*alpha/A_cap_act;
sen_index(3,:)
%SENSITIVITY WITH RESPECT TO GAMMA
gamma_sen=@(t,x_gamma)[r*x_gamma(1)-2*r*V_cap_act*x_gamma(1)-alpha*(2*V_cap_act*A_cap_act*x_gamma(1)+V_cap_act^2*x_gamma(3))-beta*(x_gamma(2)*V_cap_act+I_cap_act*x_gamma(1));theta*beta*(x_gamma(2)*V_cap_act+I_cap_act*x_gamma(1))-I_cap_act^2-2*gamma*I_cap_act*x_gamma(2)-delta*(x_gamma(2)*A_cap_act+I_cap_act*x_gamma(3));x_gamma(3)-2*A_cap_act*x_gamma(3)+x_gamma(2)*A_cap_act+I_cap_act*x_gamma(3)];
[t,x_gamma]=ode45(gamma_sen,[0 20],[0 0 0])
x_gamma_index=gamma*x_gamma
x_gamma_log_index(:,1)=x_gamma_index(:,1)/V_cap_act
x_gamma_log_index(:,2)=x_gamma_index(:,2)/I_cap_act
x_gamma_log_index(:,3)=x_gamma_index(:,3)/A_cap_act
% figure(10)
% plot(t,x_gamma_index(:,1))
% hold on
% figure(11)
% plot(t,x_gamma_index(:,2))
% hold on
% figure(12)
% plot(t,x_gamma_index(:,3))
% hold on
% figure(14)
% plot(t,x_gamma_log_index(:,1))
% hold on
% figure(15)
% plot(t,x_gamma_log_index(:,2))
% hold on
% figure(16)
% plot(t,x_gamma_log_index(:,3))
% hold on
q_gamma=size(x_gamma)
V_gamma=x_gamma(q_gamma(1),1),I_gamma=x_gamma(q_gamma(1),2),A_gamma=x_gamma(q_gamma(1),3)
sen_index(4,1)=V_gamma*gamma/V_cap_act;
sen_index(4,2)=I_gamma*gamma/I_cap_act;
sen_index(4,3)=A_gamma*gamma/A_cap_act;
sen_index(4,:)
%SENSITIVITY WITH RESPECT TO DELTA
delta_sen=@(t,x_delta)[r*x_delta(1)-2*r*V_cap_act*x_delta(1)-alpha*(2*V_cap_act*A_cap_act*x_delta(1)+V_cap_act^2*x_delta(3))-beta*(x_delta(2)*V_cap_act+I_cap_act*x_delta(1));theta*beta*(x_delta(2)*V_cap_act+I_cap_act*x_delta(1))-2*gamma*I_cap_act*x_delta(2)-I_cap_act*A_cap_act-delta*(x_delta(2)*A_cap_act+I_cap_act*x_delta(3));x_delta(3)-2*A_cap_act*x_delta(3)+x_delta(2)*A_cap_act+I_cap_act*x_delta(3)];
[t,x_delta]=ode45(delta_sen,[0 20],[0 0 0])
x_delta_index=delta*x_delta
x_delta_log_index(:,1)=x_delta_index(:,1)/V_cap_act
x_delta_log_index(:,2)=x_delta_index(:,2)/I_cap_act
x_delta_log_index(:,3)=x_delta_index(:,3)/A_cap_act
% figure(10)
% plot(t,x_delta_index(:,1))
% hold on
% figure(11)
% plot(t,x_delta_index(:,2))
% hold on
% figure(12)
% plot(t,x_delta_index(:,3))
% hold on
% figure(14)
% plot(t,x_delta_log_index(:,1))
% hold on
% figure(15)
% plot(t,x_delta_log_index(:,2))
% hold on
% figure(16)
% plot(t,x_delta_log_index(:,3))
% hold on
q_delta=size(x_delta)
V_delta=x_delta(q_delta(1),1),I_delta=x_delta(q_delta(1),2),A_delta=x_delta(q_delta(1),3)
sen_index(5,1)=V_delta*delta/V_cap_act;
sen_index(5,2)=I_delta*delta/I_cap_act;
sen_index(5,3)=A_delta*delta/A_cap_act;
sen_index(5,:)

bar(sen_index(:,2))
