clc;
clear;
close all;
format long

%parameter values that will not vary
s=1; K=1; beta=1; L=1;
%parameter values that will vary
r=.37; alpha=1e-2;
theta=.518; gamma=0.5; delta=0.08;
nu=0.05;

%main model solving in RK4 method
riparian_dim=@(t,x)[r*x(1)*(1 - x(1)/K) - alpha*(x(1).^2)*x(3) - beta*x(2)*x(1); theta*beta*x(2)*x(1) - gamma*x(2).^2 - delta*x(2)*x(3); s*x(3)*(1 - x(3)/L) + nu*x(2)*x(3)];
[t,x]=ode45(riparian_dim,[0 1000],[10 20 50]);

p_riparian=size(x);

V_cap_rev=x(p_riparian(1),1);
I_cap_rev=x(p_riparian(1),2);
A_cap_rev=x(p_riparian(1),3);

%Equilibrium Existence
a_cap = (alpha*L*nu*theta*beta)/s*(gamma+(delta*L*nu)/s);
b_cap = (r/K) + alpha*L -(alpha*L^2*nu*delta/s*(gamma + delta*L*nu/s)) + (theta*beta^2/(gamma + delta*L*nu/s));
c_cap = -r - delta*beta*L/(gamma + delta*L*nu/s);
D_cap = b_cap^2-4*a_cap*c_cap;

V1_cap = (-b_cap+sqrt(D_cap))/(2*a_cap);
V2_cap = (-b_cap-sqrt(D_cap))/(2*a_cap);

V_cap_act = V1_cap;
I_cap_act = (theta*beta*V_cap_act - delta*L) / (gamma + delta*L*nu/s);
A_cap_act = (s + nu*I_cap_act)*(L/s);

if (theta*beta*V_cap_act > delta*L)
    disp("existence of interior equilibrium satisfied")
end

%Routh-Hurwitz Stability Analysis
b1 = alpha*V_cap_act*A_cap_act + gamma*I_cap_act + s*A_cap_act/L + r*V_cap_act/K;
b2 = (r*V_cap_act*gamma*I_cap_act)/K + (r*V_cap_act*s*A_cap_act)/L + (alpha*V_cap_act*A_cap_act^2*s)/L + alpha*V_cap_act*A_cap_act*gamma*I_cap_act + (r*I_cap_act*s*A_cap_act)/L + delta*I_cap_act*nu*A_cap_act + theta*beta^2*V_cap_act*I_cap_act;
b3 = (r*V_cap_act*gamma*I_cap_act*s*A_cap_act)/K*L + (r*V_cap_act*delta*I_cap_act*nu*A_cap_act)/K + (alpha*V_cap_act*gamma*I_cap_act*s*A_cap_act^2)/L + (alpha*V_cap_act*delta*I_cap_act*nu*A_cap_act^2) + (theta*beta^2*I_cap_act*s*A_cap_act*V_cap_act)/L + alpha*V_cap_act^2*theta*beta*I_cap_act*gamma*A_cap_act;
hopf_riparian = b1*b2-b3;

if beta>=nu
    disp('condition for b1*b2-b3>0')
end

%Jacobian
a11 = r - 2*r*V_cap_act/K - 2*alpha*V_cap_act*A_cap_act - beta*I_cap_act;
a12 = -beta*V_cap_act;
a13 = -alpha*V_cap_act^2;
a21 = theta*beta*I_cap_act;
a22 = theta*beta*V_cap_act - 2*gamma*I_cap_act - delta*A_cap_act;
a23 = -delta*I_cap_act;
a31 = 0;
a32 = nu*A_cap_act;
a33= s - 2*s*A_cap_act/L + nu*I_cap_act;

J_riparian = [a11 a12 a13;a21 a22 a23;a31 a32 a33];
eigs(J_riparian)

%Figures Plotted
figure(1)
plot(t,x(:,1),'b')
hold on
figure(2)
plot(t,x(:,2),'b')
hold on
figure(3)
plot(t,x(:,3),'b')
hold on
figure(4)
plot3(x(:,3),x(:,2),x(:,1),'b')
hold on
