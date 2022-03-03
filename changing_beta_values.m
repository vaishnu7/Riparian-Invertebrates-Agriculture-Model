% clc;
% clear;
% close all;
format long
%%parameter values that will not vary
s=1; K=1; nu=1; L=1;
%%PARAMETERS WITH PERIODIC OSCILLATIONS
beta=1;
% beta=1.5;
% beta=5.5;
% beta=10.5;

r=0.9; alpha=.005;
theta=.9; gamma=.0005; delta=.003;

%main model solving in RK4 method
riparian_dim=@(t,x)[r*x(1)*(1 - (x(1)/K)) - alpha*(x(1))*x(3) - beta*x(2)*x(1); theta*beta*x(2)*x(1) - gamma*x(2)^2 - delta*x(2)*x(3); s*x(3)*(1 - (x(3)/L)) + nu*x(2)*x(3)];
[t,x]=ode45(riparian_dim,[0 1000],[.2 .3 .5]);
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
if delta<2*theta*beta && alpha<theta*beta^2
    disp("existence of interior equilibrium satisfied")
end
% if (theta*beta*V_cap_act > delta*L)
%     disp("existence of interior equilibrium satisfied")
% end
% %%dimensionless b1,b2,b3
% b1_dim = alpha*V_cap_act*A_cap_act + gamma*I_cap_act + A_cap_act + r*V_cap_act;
% b2_dim = r*V_cap_act*gamma*I_cap_act + r*V_cap_act*A_cap_act + alpha*V_cap_act*(A_cap_act^2) + (alpha*V_cap_act*A_cap_act*gamma*I_cap_act) + gamma*I_cap_act*A_cap_act + delta*I_cap_act*nu*A_cap_act + theta*V_cap_act*I_cap_act;
% b3_dim = r*V_cap_act*gamma*I_cap_act*A_cap_act + r*V_cap_act*delta*I_cap_act*nu*A_cap_act + alpha*V_cap_act*gamma*I_cap_act*(A_cap_act^2) + alpha*V_cap_act*delta*I_cap_act*nu*A_cap_act^2 + theta*I_cap_act*A_cap_act*V_cap_act + alpha*V_cap_act^2*theta*I_cap_act*nu*A_cap_act;
% hopf_riparian_dim = b1_dim*b2_dim-b3_dim
% 
% %Routh-Hurwitz Stability Analysis
b1 = alpha*V_cap_act*A_cap_act + gamma*I_cap_act + (s*A_cap_act/L) + (r*V_cap_act/K)
b2 = ((r*V_cap_act*gamma*I_cap_act)/K) + ((r*V_cap_act*s*A_cap_act)/L) + ((alpha*V_cap_act*(A_cap_act^2)*s)/L) + (alpha*V_cap_act*A_cap_act*gamma*I_cap_act) + ((gamma*I_cap_act*s*A_cap_act)/L) + (delta*I_cap_act*nu*A_cap_act) + (theta*beta^2*V_cap_act*I_cap_act);
b3 = ((r*V_cap_act*gamma*I_cap_act*s*A_cap_act)/(K*L)) + ((r*V_cap_act*delta*I_cap_act*nu*A_cap_act)/K) + ((alpha*V_cap_act*gamma*I_cap_act*s*A_cap_act^2)/L) + (alpha*V_cap_act*delta*I_cap_act*nu*A_cap_act^2) + ((theta*beta^2*I_cap_act*s*A_cap_act*V_cap_act)/L) + alpha*V_cap_act^2*theta*beta*I_cap_act*nu*A_cap_act;
hopf_riparian = b1*b2-b3
% 
if hopf_riparian>0
    disp('condition for b1*b2-b3>0')
end
% %Jacobian
% a11_1= - (r*V_cap_act/K)- alpha*V_cap_act*A_cap_act
a11 = r - (2*r*V_cap_act/K) - 2*alpha*V_cap_act*A_cap_act - beta*I_cap_act;
a12 = -beta*V_cap_act;
a13 = -alpha*V_cap_act^2;
a21 = theta*beta*I_cap_act;
a22 = theta*beta*V_cap_act - 2*gamma*I_cap_act - delta*A_cap_act;
% a22=-gamma*I_cap_act
a23 = -delta*I_cap_act;
a31 = 0;
a32 = nu*A_cap_act;
% a33=-s*A_cap_act/L
a33= s - (2*s*A_cap_act/L) + nu*I_cap_act;
J_riparian = [a11 a12 a13;a21 a22 a23;a31 a32 a33];
eigs(J_riparian)
% %Global Stability Condition (using Poincare Bendixosn theorem)
if gamma>delta && r>alpha*K^2
    disp('condition 1 for Global Stability satisfied')
end
if beta<r/K
    disp('condition 2 for Global Stability satisfied')
end
% if nu<s/L && beta<r/K
%     disp('condition 2 for Global Stability satisfied')
% end
if theta*beta<gamma
    disp('condition 3 for Global Stability satisfied')
end
% %Global Stability Condition (using Lyapunov)
if alpha^2*K^2<4*(delta*s*r)/(L*theta*nu*K)
    disp('condition 4 (Lyapunov) for Global Stability satisfied')
end
%%tranvsersality condition with respect to delta
trans_value=theta*I_cap_act*V_cap_act*(alpha*V_cap_act*A_cap_act-2*gamma*beta*I_cap_act-2*r*beta*V_cap_act-2*alpha*beta*V_cap_act*A_cap_act)
if trans_value~= 0
      disp('transversality condition satisfied')
end 

%Figures Plotted
figure(1)
plot(t,x(:,1),'b')
xlabel('time') 
ylabel('riparian vegetation')
label1 = '$\beta=0.05$';
label2 = '$\beta=1.5$';
label3 = '$\beta=5.5$';
legend(label1,label2,label3,'Interpreter','latex')
hold on
figure(2)
plot(t,x(:,2),'b')
xlabel('time') 
ylabel('terrestrial invertebrates')
label1 = '$\beta=0.05$';
label2 = '$\beta=1.5$';
label3 = '$\beta=5.5$';
legend(label1,label2,label3,'Interpreter','latex')
hold on
figure(3)
plot(t,x(:,3),'b')
xlabel('time') 
ylabel('agriculture')
label1 = '$\beta=0.05$';
label2 = '$\beta=1.5$';
label3 = '$\beta=5.5$';
legend(label1,label2,label3,'Interpreter','latex')
hold on
figure(4)
plot3(x(:,1),x(:,2),x(:,3),'b')
xlabel('riparian vegetation') 
ylabel('terrestrial invertebrates')
zlabel('agriculture')
label1 = '$\beta=0.05$';
label2 = '$\beta=1.5$';
label3 = '$\beta=5.5$';
legend(label1,label2,label3,'Interpreter','latex')
hold on
% figure(5)
% plot(x(:,2),x(:,1),'k')
% xlabel('terrestrial invertebrates') 
% ylabel('riparian vegetation') 
% hold on