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
delT1=(T12max-T12min)/50.;
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

%%
var1=zeros(nruns_tot,1);
var2=zeros(nruns_tot,1);
var3=zeros(nruns_tot,1);
var4=zeros(nruns_tot,1);
var5=zeros(nruns_tot,1);
var6=zeros(nruns_tot,1);
var7=zeros(nruns_tot,1);
var8=zeros(nruns_tot,1);
%Max_evalue_ratio=zeros(nruns_tot,1);

%outside loop over T1 values of 2nd component; inside loop over D values of 2nd component
SigmaMat2=zeros(nrunsT1,nrunsD,8);
for kk=1:nrunsT1
    x7=x7v(kk);
    parfor jj=1:nrunsD
        x8=x8v(jj);
        P_3D=[x1 x2 x3 x4 x5 x6 x7 x8];  %Vector of unknown parameters
        A=Derivatives_3D(TE,tau,b,P_3D);
        HLM=A'*A;
        Cov=inv(HLM);
        Sigma=diag(Cov);
        SigmaMat2(kk,jj,:)=Sigma*STDnoise^2;

    end
end
%%
[X,Y] = meshgrid(x8v, x7v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
surf(X,Y,sqrt(SigmaMat2(:,:,1)),'FaceAlpha',0.8)
ylabel('  T_{1,2} ','FontWeight','bold','FontSize',44)
xlabel('ADC_{2} ','FontWeight','bold','FontSize',44)
zlabel({'Std Dev of Derived'; 'c_1 Value'},'FontWeight','bold','FontSize',44)
axis([ADC2min ADC2max T12min T12max]) 
% zlim([0 1700])
f = gcf;
f.Position = [100 100 700 500];

h2 = figure;
surf(X,Y,sqrt(SigmaMat2(:,:,5)),'FaceAlpha',0.8)
ylabel('  T_{1,2} ','FontWeight','bold','FontSize',44)
xlabel('ADC_{2} ','FontWeight','bold','FontSize',44)
zlabel({'Std Dev of Derived'; 'c_2 Value'},'FontWeight','bold','FontSize',44)
axis([ADC2min ADC2max T12min T12max]) 
% zlim([0 1700])
f = gcf;
f.Position = [100 100 700 500];

h3 = figure;
surf(X,Y,sqrt(SigmaMat2(:,:,2)),'FaceAlpha',0.8)
ylabel('T_{1,2} ','FontWeight','bold','FontSize',44)
xlabel('ADC_{2} ','FontWeight','bold','FontSize',44)
zlabel({'Std Dev of Derived'; 'T_{2,1} Value'},'FontWeight','bold','FontSize',44)
axis([ADC2min ADC2max T12min T12max]) 
% zlim([0 1700])
f = gcf;
f.Position = [100 100 700 500];

h4 = figure;
surf(X,Y,sqrt(SigmaMat2(:,:,6)),'FaceAlpha',0.8)
ylabel('  T_{1,2} ','FontWeight','bold','FontSize',44)
xlabel('ADC_{2} ','FontWeight','bold','FontSize',44)
zlabel({'Std Dev of Derived'; 'T_{2,2} Value'},'FontWeight','bold','FontSize',44)
axis([ADC2min ADC2max T12min T12max]) 
% zlim([0 1700])

f = gcf;
f.Position = [100 100 700 500];
%% save fig

saveas(h1,'analytic_3D_c1_std','epsc')
saveas(h2,'analytic_3D_c2_std','epsc')
saveas(h3,'analytic_3D_T2A_std','epsc')
saveas(h4,'analytic_3D_T2B_std','epsc')


