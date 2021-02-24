%% 1D analytic
clc;clear;close all
%% Parameter settings
%%c1=offset term; c1, c2=populations; tau1, tau2=time constants;
%%TE=echo times, which are the sampling times
%%model:  M(x,TE)=c1*exp(-TE/tau1)+c2*exp(-TE/tau2)
c1=0.3;
c2=(1-c1);
tau1=60;

tau2min=30;
tau2max=100;
% del=(tau2max-tau2min)/100.;
del=0.1;
tau2v=tau2min:del:tau2max;
%holesizeleft is the 1/2-range of values to the left of tau1 that are omitted 
%holesizeright is the 1/2-range of values to the right of tau1 that are
%omitted
holesizeleft=0;
holesizeright=0;

i1=find(tau2v>(tau1-holesizeleft),1,'first');
i2=find(tau2v>(tau1+holesizeright),1,'first');

nruns=length(tau2v);

TE=(8:8:512)';
%%
%defining the number of sub-plots

plotrows=2;
plotcols=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var1=zeros(nruns,1);
var2=zeros(nruns,1);
var3=zeros(nruns,1);
var4=zeros(nruns,1);
CJ=zeros(nruns,1);
CJTJ=zeros(nruns,1);
SNR = 800;
STDnoise=(1/SNR);
%%model:  M(x,TE)=c1*exp(-TE/tau1)+c2*exp(-TE/tau2)
for kk=1:nruns
    tau2=tau2v(kk);
    P_biexp=[c1 c2 tau1 tau2];  %Vector of unknown parameters

%%%Derivative matrix computation

    J=-T2Derivatives_biexp_no_const(TE,P_biexp);%keep this line as is, including "-", for non-regularized. It works. 
%**********************
%%%%(see pg 155, Hansen Least Squares).  (This was MMB's Hessian.)
    Hessian_approx=J'*J;
%%%Covariance matrix computation, assuming unit noise standard deviation
    Cov=Hessian_approx^(-1);
    Sigmas=diag(Cov);
%%model:  M(x,TE)=c1+c2*exp(-TE/tau1)+c3*exp(-TE/tau2)
%var3 and var4 are the variances of the two T2 values--see def of P_biexp
    var3(kk)=Sigmas(3)*STDnoise^2;
    var4(kk)=Sigmas(4)*STDnoise^2;
%var1 and var2 are the variances of the two population fractions--see definition of
%P_biexp
    var1(kk)=Sigmas(1)*STDnoise^2;
    var2(kk)=Sigmas(2)*STDnoise^2;
%   Max_evalue_ratio(kk)=max(Norm_Eig_Hessian_approx)/min(Norm_Eig_Hessian_approx);
CJ(kk)=cond(J);
CJTJ(kk)=cond(Hessian_approx);
end

%% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = figure;

%plotting sqrt(var)=SD of parameters

%%% T2A
% subplot(plotrows,plotcols,1) 
%plot(tau2v,var4,'LineWidth',3)
% var3=var3/(min(vertcat(var3(1:i1), var3(i2:length(var3)))));
semilogy(tau2v(1:i1),sqrt(var3(1:i1)),'-k','LineWidth',3);
hold on
semilogy(tau2v(i2:length(tau2v)),sqrt(var3(i2:length(var3))),'-k','LineWidth',3)


xlabel('T_{2,2}','FontWeight','bold','FontSize',40)

ylabel({'Std Dev of Derived'; 'T_{2,1} Value'},'FontWeight','bold','FontSize',40)
ylim([1e0 1e5])
hold off
saveas(ax1,'analytic_1D_T2A_std','epsc')

%%
%%% T2B
ax2 = figure;
% subplot(plotrows,plotcols,2)
% var4=var4/(min(vertcat(var4(1:i1), var4(i2:length(var4)))));
semilogy(tau2v(1:i1),sqrt(var4(1:i1)),'-k','LineWidth',3);
hold on
semilogy(tau2v(i2:length(tau2v)),sqrt(var4(i2:length(var3))),'-k','LineWidth',3)
xlabel('T_{2,2}','FontWeight','bold','FontSize',40)

ylabel({'Std Dev of Derived'; 'T_{2,2} Value'},'FontWeight','bold','FontSize',40)
ylim([1e0 1e5])
saveas(ax2,'analytic_1D_T2B_std','epsc')

%%% c1
ax3 = figure;
% subplot(plotrows,plotcols,3)
% var1=var1/(min(vertcat(var1(1:i1), var1(i2:length(var4)))));
semilogy(tau2v(1:i1),sqrt(var1(1:i1)),'-k','LineWidth',3);
hold on
semilogy(tau2v(i2:length(tau2v)),sqrt(var1(i2:length(var3))),'-k','LineWidth',3)
xlabel('T_{2,2}','FontWeight','bold','FontSize',40)

ylabel({'Std Dev of Derived'; 'c_1 Value'},'FontWeight','bold','FontSize',40)
ylim([1e-2 1e6])
saveas(ax3,'analytic_1D_c1_std','epsc')

%%% c2
ax4 = figure;
% subplot(plotrows,plotcols,4)
% var2=var2/(min(vertcat(var2(1:i1), var2(i2:length(var4)))));
semilogy(tau2v(1:i1),sqrt(var2(1:i1)),'-k','LineWidth',3);
hold on
semilogy(tau2v(i2:length(tau2v)),sqrt(var2(i2:length(var3))),'-k','LineWidth',3)
xlabel('T_{2,2}','FontWeight','bold','FontSize',40)

ylabel({'Std Dev of Derived'; 'c_2 Value'},'FontWeight','bold','FontSize',40)

ylim([1e-2 1e6])

% f = gcf;
% f.Position = [100 100 1500 900];
hold off
saveas(ax4,'analytic_1D_c2_std','epsc')

%%
holesizeleft=0.1;
holesizeright=0;
% holesizeleft=22;
% holesizeright=200;
i11=find(tau2v>(tau1-holesizeleft),1,'first');
i21=find(tau2v>(tau1+holesizeright),1,'first');

ax5 = figure;
% subplot(plotrows,plotcols,3)
semilogy(tau2v(1:i11),CJ(1:i11),'-k','LineWidth',3)
hold on
semilogy(tau2v(i21:length(tau2v)),CJ(i21:length(CJ)),'-k','LineWidth',3)
% xlabel('T_{2,B}','FontSize',18)
% ylabel({'Condition number of'; 'Jacobian matrix'},'FontSize',18)


xlabel('T_{2,2}','FontWeight','bold','FontSize',30)

ylabel({'Condition Number of'; 'Jacobian Matrix'},'FontWeight','bold','FontSize',30)

hold off
saveas(ax5,'analytic_1D_jacobian_cond','epsc')

