%% 2D, analytic,
% Apr 1 2020
clear;clc;close all

%% Input parameters
%%TE=echo times, which are the direct dimension sampling times
%%tau=indirect dimension sampling times
%%model:  s(x,t_i, tau_j)=c1*exp(-t_i/T2A)*(1-2exp(-tau_j/xT1A))+c2*exp(-t_i/T2B)*(1-2exp(-tau_j/T1B))
c1=0.3; c2=(1-c1);
T2A=60; T2B=45;
T1A=1000;

T1Bmin=900;
T1Bmax=1100;
del=1;
T1Bv=T1Bmin:del:T1Bmax; % This is T1B values
nruns=length(T1Bv);
TE=(8:8:512)';
tau=50:200:5000; % This is TI
SNR = 800;
STDnoise=(1/SNR);
%%
var1=zeros(nruns,1);
var2=zeros(nruns,1);
var3=zeros(nruns,1);
var4=zeros(nruns,1);
var5=zeros(nruns,1);
var6=zeros(nruns,1);

CA=zeros(nruns,1);
% Max_evalue_ratio=zeros(nruns,1);
aa = [];
for kk=1:nruns

    T1B=T1Bv(kk);
    %     T1A=T1B;
    P_2D=[c1 T2A T1A c2 T2B T1B]; 

%%
%%%Derivative matrix computation
%%A is sort of a Jacobian matrix
    A=T1T2Derivatives_2D(TE,tau,P_2D);
%HLM means "Hessian-like matrix", that is, A'*A is Hessian-like, but
%not a Hessian (as far as I can tell)
    HLM=A'*A;
%%%(see pg 155). 
%     Norm_Eig_HLM=eig(HLM)/max(eig(HLM));
%%%Covariance matrix computation, assuming unit noise standard deviation
    Cov=inv(HLM);
    Sigmas=diag(Cov);
%%model:  s(x,t_i, tau_j)=x1*exp(-t_i/x2)*(1-2exp(-tau_j/x3))+x4*exp(-t_i/x5)*(1-2exp(-tau_j/x6))
%var2 and var5 are the variances of the two T2 values
    var2(kk)=Sigmas(2)*STDnoise^2;
    var5(kk)=Sigmas(5)*STDnoise^2;
%var1 and var4 are the variances of the two component size values
    var1(kk)=Sigmas(1)*STDnoise^2;
    var4(kk)=Sigmas(4)*STDnoise^2;
%     Max_evalue_ratio(kk)=max(Norm_Eig_HLM)/min(Norm_Eig_HLM);
%var3 and var6 are the variances of the two T1 values
    var3(kk)=Sigmas(3)*STDnoise^2;
    var6(kk)=Sigmas(6)*STDnoise^2;
    aa = [aa;Sigmas(4)/Sigmas(5)];
    CA(kk)=cond(A);
end

%% plots
%var1=c1,var2=T2A,var3=T1A
%var4=c2,var5=T2B,var6=T1B

ax1 = figure;

plot(T1Bv,sqrt(var2),'LineWidth',3,'color','k')
% semilogy(x6v,sqrt(var2),'LineWidth',3)
xlabel('T_{1,2}','FontWeight','bold','FontSize',48)
ylabel({'Std Dev of Derived'; 'T_{2,1} Value'},'FontWeight','bold','FontSize',48)
% ylim([0 2.5])

f = gcf;
f.Position = [100 100 700 500];
saveas(ax1,'analytic_2D_T2A_std','epsc')

ax2 = figure;
plot(T1Bv,sqrt(var5),'LineWidth',3,'color','k')
xlabel('T_{1,2}','FontWeight','bold','FontSize',48)
ylabel({'Std Dev of Derived'; 'T_{2,2} Value'},'FontWeight','bold','FontSize',48)
% ylim([0 2.5])

f = gcf;
f.Position = [100 100 700 500];
saveas(ax2,'analytic_2D_T2B_std','epsc')

%%

% figure
% plot(T1Bv,sqrt(var3),'LineWidth',3,'color','k')
% xlabel('T_{1,2}','FontWeight','bold','FontSize',30)
% ylabel({'Std Dev of Derived'; 'T_{1,1} Value'},'FontWeight','bold','FontSize',30)
% f = gcf;
% f.Position = [100 100 700 500];
% 
% figure;
% plot(T1Bv,sqrt(var6),'LineWidth',3,'color','k')
% xlabel('T_{1,2}','FontWeight','bold','FontSize',30)
% ylabel({'Std Dev of Derived'; 'T_{1,2} Value'},'FontWeight','bold','FontSize',30)
% f = gcf;
% f.Position = [100 100 700 500];
%%

ax3 = figure;
plot(T1Bv,sqrt(var1),'LineWidth',3,'color','k')
% semilogy(x6v,sqrt(var2),'LineWidth',3)
xlabel('T_{1,2}','FontWeight','bold','FontSize',48)
ylabel({'Std Dev of Derived'; 'c_1 Value'},'FontWeight','bold','FontSize',48)
% ylim([0 0.2])

f = gcf;
f.Position = [100 100 700 500];
saveas(ax3,'analytic_2D_c1_std','epsc')

ax4 = figure;
plot(T1Bv,sqrt(var4),'LineWidth',3,'color','k')
xlabel('T_{1,2}','FontWeight','bold','FontSize',48)
ylabel({'Std Dev of Derived'; 'c_2 Value'},'FontWeight','bold','FontSize',48)
% ylim([0 0.2])
f = gcf;
f.Position = [100 100 700 500];
saveas(ax4,'analytic_2D_c2_std','epsc')

%%
% 
ax5 = figure;
plot(T1Bv,CA,'LineWidth',3,'color','k')
xlabel('T_{1,2}','FontWeight','bold','FontSize',36)
ylabel({'Condition Number of', 'Jacobian Matrix'},'FontWeight','bold','FontSize',36)

f = gcf;
f.Position = [100 100 700 500];

saveas(ax5,'analytic_2D_jacobian','epsc')
