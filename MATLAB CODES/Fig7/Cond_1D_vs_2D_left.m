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

T1Bmin=800;
T1Bmax=1200;
del=1;
T1Bv=T1Bmin:del:T1Bmax; % This is T1B values
nruns=length(T1Bv);
TE=(8:8:512)';
tau=1:200:5000; % This is TI
%%
var1=zeros(nruns,1);
var2=zeros(nruns,1);
var3=zeros(nruns,1);
var4=zeros(nruns,1);
var5=zeros(nruns,1);
var6=zeros(nruns,1);

CA=zeros(nruns,1);
% Max_evalue_ratio=zeros(nruns,1);
for kk=1:nruns

    T1B=T1Bv(kk);
    if T1B == T1A
        P_2D=[c1 T2A T1A c2 T2B T1B]; 
        AA=T1T2Derivatives_2D(TE,tau,P_2D);
    end
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
    Cov=pinv(HLM);
    Sigmas=diag(Cov);
%%model:  s(x,t_i, tau_j)=x1*exp(-t_i/x2)*(1-2exp(-tau_j/x3))+x4*exp(-t_i/x5)*(1-2exp(-tau_j/x6))
%var2 and var5 are the variances of the two T2 values
    var2(kk)=Sigmas(2);
    var5(kk)=Sigmas(5);
%var1 and var4 are the variances of the two component size values
    var1(kk)=Sigmas(1);
    var4(kk)=Sigmas(4);
%     Max_evalue_ratio(kk)=max(Norm_Eig_HLM)/min(Norm_Eig_HLM);
%var3 and var6 are the variances of the two T1 values
    var3(kk)=Sigmas(3);
    var6(kk)=Sigmas(6);
    CA(kk)=cond(A);
end

%% 1D condition number
P_biexp=[c1 c2 T2A T2B];  %Vector of unknown parameters
J=-T2Derivatives_biexp_no_const(TE,P_biexp);
CJ = cond(J);
CJ = ones(size(T1Bv,2),1) * CJ;
%%
% 
figure
hold on 
plot(T1Bv, CJ,'LineWidth',4);
plot(T1Bv,CA,'LineWidth',4)
ylim([1 3]*1e4) % left panel [0.3, 0.7, 60, 45]
% ylim([5.6 7.6]*1e3) % right panel [0.2, 0.8, 40, 90]

legend({'1D','2D'}, 'FontSize',36)
xlabel('T_{1,2}','FontWeight','bold','FontSize',44)
ylabel({'Condition Number of'; 'Jacobian Matrix'},'FontWeight','bold','FontSize',44)
% title(['c1 = ',num2str(c1), ', c2 = ',num2str(c2),', T21 = ',num2str(T2A), ', T22 = ', num2str(T2B), ', T11 = ', num2str(T1A), ', T12 varies'], 'FontSize',18)

f = gcf;
f.Position = [100 100 700 500];

%% save fig
saveas(f,'analytic_1D_vs_2D_jacobian_cond','epsc')