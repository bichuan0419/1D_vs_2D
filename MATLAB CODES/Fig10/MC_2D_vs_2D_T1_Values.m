%1D vs 2D
% This is by fixing T1A = 3*T1B
clear;clc;
close all;
c1=0.3;c2=0.7;
TE=8:8:512;
TI = 50:200:5000;
nruns=1000;

opt = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt');
T21=45;T22=60;

nbins=20;

Z=c1*exp(-TE/T21) + c2*exp(-TE/T22);
ZMAX=max(max(abs(Z)));
SNR=400;
STDnoise=ZMAX*(1/SNR);

%%%%%%%%%re-run with different values for T11 and T12
%% Now run the 2D model for comparison with the above 1D model
[X,Y] = meshgrid(TE,TI);
XY(:,:,1) = X;
XY(:,:,2) = Y;
% % Create Objective Function; this is the fitting model: 
surfit = @(B,XY) B(5)*exp(-TE/B(1)) .* (1-2.*exp((-TI/B(2))))'+ B(6)*exp(-TE/B(3)) .* (1-2.*exp((-TI/B(4))))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment out the below since want to do the comparison based on the same
% noise (with the sqrt(length(TI)) term still required for equal time)
% ZMAX=max(max(abs(Z)));
% STDnoise=ZMAX*(1/SNR)

%% 2D with unequal T1 values: T1A = 1000, T1B = 200
T11=300;T12=900;
Z=c1*exp(-TE/T21).*(1-2.*exp(-TI/T11))' + c2*exp(-TE/T22).*(1-2.*exp(-TI/T12))';

TT21_2DA = zeros(nruns,1);
TT22_2DA = zeros(nruns,1);
TT11_2DA = zeros(nruns,1);
TT12_2DA = zeros(nruns,1);
c1_2DA = zeros(nruns,1);
c2_2DA = zeros(nruns,1);

parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T11+100*rand() T22+10*rand()...
    T12+100*rand() rand() rand()], XY, ZN, [],[],opt);
TT21_2DA(kk)=BN(1);TT11_2DA(kk)=BN(2);TT22_2DA(kk)=BN(3);TT12_2DA(kk)=BN(4);
c1_2DA(kk) = BN(5); c2_2DA(kk) = BN(6);
end
for i = 1:nruns
    if TT21_2DA(i) > TT22_2DA(i)
        a = TT21_2DA(i);
        b = TT22_2DA(i);
        TT21_2DA(i) = b;
        TT22_2DA(i) = a;
        
        c = c1_2DA(i);
        d = c2_2DA(i);
        c1_2DA(i) = d;
        c2_2DA(i) = c;
    end
end
%%% plot
h1 = figure;
tiledlayout(1,2)

ax3 = nexttile;
histogram(TT21_2DA, nbins)
hold on
histogram(TT22_2DA, nbins)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_2DA),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_2DA),2)]},'FontSize',30, 'Location', 'North') 
% title( {['2D: T_{1,1}=' num2str(T11), ', T_{1,2}=' num2str(T12)]},'FontSize',30); 

ax4 = nexttile;

histogram(c1_2DA, nbins)
hold on
histogram(c2_2DA, nbins)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold') 
legend({['c_1', ', SD = ', num2str(std(c1_2DA),2)],['c_2', ', SD = ', num2str(std(c2_2DA),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: T_{1,1} = ' num2str(T11), ', T_{1,2} = ' num2str(T12)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1845         159        1284         348])


%% 2D with equal T1 values: T1A = 1000, T1B = 1200
T11=600;T12=1800;
Z=c1*exp(-TE/T21).*(1-2.*exp(-TI/T11))' + c2*exp(-TE/T22).*(1-2.*exp(-TI/T12))';

TT21_2DB = zeros(nruns,1);
TT22_2DB = zeros(nruns,1);
TT11_2DB = zeros(nruns,1);
TT12_2DB = zeros(nruns,1);
c1_2DB = zeros(nruns,1);
c2_2DB = zeros(nruns,1);

parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T11+100*rand() T22+10*rand()...
    T12+100*rand() rand() rand()], XY, ZN, [],[],opt);
TT21_2DB(kk)=BN(1);TT11_2DB(kk)=BN(2);TT22_2DB(kk)=BN(3);TT12_2DB(kk)=BN(4);
c1_2DB(kk) = BN(5); c2_2DB(kk) = BN(6);
end
for i = 1:nruns
    if TT21_2DB(i) > TT22_2DB(i)
        a = TT21_2DB(i);
        b = TT22_2DB(i);
        TT21_2DB(i) = b;
        TT22_2DB(i) = a;
        
        c = c1_2DB(i);
        d = c2_2DB(i);
        c1_2DB(i) = d;
        c2_2DB(i) = c;
    end
end
%%% plot
h2 = figure;
tiledlayout(1,2)

ax5 = nexttile;
histogram(TT21_2DB, nbins)
hold on
histogram(TT22_2DB, nbins)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_2DB),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_2DB),2)]},'FontSize',30, 'Location', 'North') 
% title( {['2D: T_{1,1}=' num2str(T11), ', T_{1,2}=' num2str(T12)]},'FontSize',30); 

ax6 = nexttile;

histogram(c1_2DB, nbins)
hold on
histogram(c2_2DB, nbins)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold') 
legend({['c_1', ', SD = ', num2str(std(c1_2DB),2)],['c_2', ', SD = ', num2str(std(c2_2DB),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: T_{1,1} = ' num2str(T11), ', T_{1,2} = ' num2str(T12)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1845         159        1284         348])


%% 2D with unequal T1 values: T1A = 1000, T1B = 2000
T11=900;T12=2700;
Z=c1*exp(-TE/T21).*(1-2.*exp(-TI/T11))' + c2*exp(-TE/T22).*(1-2.*exp(-TI/T12))';

TT21_2DC = zeros(nruns,1);
TT22_2DC = zeros(nruns,1);
TT11_2DC = zeros(nruns,1);
TT12_2DC = zeros(nruns,1);
c1_2DC = zeros(nruns,1);
c2_2DC = zeros(nruns,1);

parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T11+100*rand() T22+10*rand()...
    T12+100*rand() rand() rand()], XY, ZN, [],[],opt);
TT21_2DC(kk)=BN(1);TT11_2DC(kk)=BN(2);TT22_2DC(kk)=BN(3);TT12_2DC(kk)=BN(4);
c1_2DC(kk) = BN(5); c2_2DC(kk) = BN(6);
end
for i = 1:nruns
    if TT21_2DC(i) > TT22_2DC(i)
        a = TT21_2DC(i);
        b = TT22_2DC(i);
        TT21_2DC(i) = b;
        TT22_2DC(i) = a;
        
        c = c1_2DC(i);
        d = c2_2DC(i);
        c1_2DC(i) = d;
        c2_2DC(i) = c;
    end
end
%%% plot
h3 = figure;
tiledlayout(1,2)

ax7 = nexttile;
histogram(TT21_2DC, nbins)
hold on
histogram(TT22_2DC, nbins)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_2DC),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_2DC),2)]},'FontSize',30, 'Location', 'North') 
% title( {['2D: T_{1,1}=' num2str(T11), ', T_{1,2}=' num2str(T12)]},'FontSize',30); 

ax8 = nexttile;

histogram(c1_2DC, nbins)
hold on
histogram(c2_2DC, nbins)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold') 
legend({['c_1', ', SD = ', num2str(std(c1_2DC),2)],['c_2', ', SD = ', num2str(std(c2_2DC),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: T_{1,1} = ' num2str(T11), ', T_{1,2} = ' num2str(T12)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1845         159        1284         348])


linkaxes([ax3,ax5,ax7],'xy')
linkaxes([ax4,ax6,ax8],'xy')
% 
%% save figures
saveas(h1,'MC_2D_vs_2D_A','epsc')
saveas(h2,'MC_2D_vs_2D_B','epsc')
saveas(h3,'MC_2D_vs_2D_C','epsc')

