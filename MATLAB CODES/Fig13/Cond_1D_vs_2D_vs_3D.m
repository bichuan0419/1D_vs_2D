clear;clc;
close all;
%%
%Input parameters
%%TE=echo times, which are the direct dimension sampling times
%%tau=indirect dimension sampling times
%%b=2nd indirect dimension sampling times
%%model: s(x,t_i,tau_j,b_k)=x1*exp(-t_i/x2)*(1-2exp(-tau_j/x3))*exp(-b_k*x4)+x5*exp(-t_i/x6)*(1-2exp(-tau_j/x7))*exp(-b_k*x8)

%component 1:  fraction, T2, T1, D
c1=0.3;
% T21=20.;
T21=45.;
T11=1000;
ADC1=1.50;


x1=c1;
x2=T21;
x3=T11;
x4=ADC1;

%component 2:  fraction, T2, T1 range, D range
c2=(1-c1);
T22=60;
T12min=800;
T12max=1200;
delT1=(T12max-T12min)/100.;
ADC2min=0.5;
ADC2max=4;
delADC2=0.1;
x5=c2;
x6=T22;
x7v=T12min:delT1:T12max;
x8v=ADC2min:delADC2:ADC2max;

nrunsT1=length(x7v);
nrunsD=length(x8v);
nruns_tot=nrunsT1*nrunsD;

%vectors of measurement "times"--echo times, inversion times, b values

TE=8:8:512;
tau=50:200:5000;
b=0:0.25:2;
nTE = length(TE);
nTI = length(tau);
nADC = length(b);

SNR = 200;
STDnoise=(1/SNR);
%% 1D condition number

P_biexp=[c1 c2 T21 T22];  %Vector of unknown parameters
J=-T2Derivatives_biexp_no_const(TE,P_biexp);
CJ = cond(J);
CJ_1D = ones(nrunsT1,nrunsD) * CJ;

%% 2D condition number
CJ_2D=zeros(nrunsT1,nrunsD);

for kk=1:nrunsT1

    T12=x7v(kk);
%     if T12 == T11
%         P_2D=[c1 T21 T11 c2 T22 T12]; 
%         AA=T1T2Derivatives_2D(TE,tau,P_2D);
%     end
    %     T1A=T1B;
    P_2D=[c1 T21 T11 c2 T22 T12]; 

%%
%%%Derivative matrix computation
%%A is sort of a Jacobian matrix
    A=T1T2Derivatives_2D(TE,tau,P_2D);

    CJ_2D(kk,:)=cond(A)*ones(1,nrunsD);
end

%% 3D condition number
%outside loop over T1 values of 2nd component; inside loop over D values of 2nd component
CJ_3D = zeros(nrunsT1,nrunsD);
for kk=1:nrunsT1
    x7=x7v(kk);
    parfor jj=1:nrunsD
        x8=x8v(jj);
        P_3D=[x1 x2 x3 x4 x5 x6 x7 x8];  %Vector of unknown parameters
        A=Derivatives_3D(TE,tau,b,P_3D);
        CJ_3D(kk,jj)= cond(A);

    end
end
%%
[X,Y] = meshgrid(x8v, x7v);
%%
h1=figure;
surf(X,Y,CJ_1D,'FaceAlpha',0.2,'FaceColor','red')
hold on
surf(X,Y,CJ_2D,'FaceAlpha',0.6,'FaceColor','green')

surf(X,Y,CJ_3D,'FaceAlpha',1,'FaceColor','blue')
legend({'1D','2D','3D'},'FontSize',24)
ylabel('T_{1,2} ','FontWeight','bold','FontSize',30)
xlabel('ADC_{2} ','FontWeight','bold','FontSize',30)
zlabel({'Condition Number of', 'Jacobian Matrix'},'FontWeight','bold','FontSize',30)
axis([ADC2min ADC2max T12min T12max]) 
set(gcf,'position',[440   168   829   630])
view([-50 15])

%% plot slice at T11 = T12
T1_equal_idx = find(x7v == T11);
from_ADC_1D = CJ_1D(T1_equal_idx,:);
from_ADC_2D = CJ_2D(T1_equal_idx,:);
from_ADC_3D = CJ_3D(T1_equal_idx,:);
h2 = figure;
hold on
plot(x8v,from_ADC_1D, 'LineWidth',6,'Color','red');
plot(x8v,from_ADC_2D, 'LineWidth',3,'Color','green');
plot(x8v,from_ADC_3D, 'LineWidth',3,'Color','blue');
xlabel('ADC_{2} ','FontWeight','bold','FontSize',36)
ylabel({'Condition Number of', 'Jacobian Matrix'},'FontWeight','bold','FontSize',36)
legend({'1D','2D','3D'},'FontSize',30)
title('Side View at T_{1,1} = T_{1,2}','FontWeight','bold','FontSize',36)
%% plot slice at ADC1 = ADC2
ADC_equal_idx = find(x8v == ADC1);
from_T1_1D = CJ_1D(:,ADC_equal_idx);
from_T1_2D = CJ_2D(:,ADC_equal_idx);
from_T1_3D = CJ_3D(:,ADC_equal_idx);
h3 = figure;
hold on 
plot(x7v,from_T1_1D, 'LineWidth',3,'Color','red');
plot(x7v,from_T1_2D, 'LineWidth',6,'Color','green');
plot(x7v,from_T1_3D, 'LineWidth',3,'Color','blue');
xlabel('T_{1,2} ','FontWeight','bold','FontSize',36)
ylim([1 3]*1e4)
ylabel({'Condition Number of', 'Jacobian Matrix'},'FontWeight','bold','FontSize',36)
legend({'1D','2D','3D'},'FontSize',30)
title('Side View at ADC_1 = ADC_2','FontWeight','bold','FontSize',36)

%% save fig
saveas(h1,'analytic_1D_vs_2D_vs_3D_cond','epsc')
saveas(h2,'analytic_1D_vs_2D_vs_3D_cond_from_ADC','epsc')
saveas(h3,'analytic_1D_vs_2D_vs_3D_cond_from_T1','epsc')