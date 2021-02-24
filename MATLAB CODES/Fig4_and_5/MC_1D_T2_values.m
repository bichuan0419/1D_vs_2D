%1D vs 2D
% This is by fixing T1A = 3*T1B
clear;clc;
close all;
c1=0.7;c2=1-c1;
TE=8:8:512;
nruns=1000;

opt = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt');

SNR=10000;
sd_c1 = zeros(3,2);
sd_c2 = zeros(3,2);
sd_T21 = zeros(3,2);
sd_T22 = zeros(3,2);

%% T21 = 60, T22 = 45

nbins_T2 = linspace(30,100,100);
nbins_c = linspace(0,1,100);

T21 = 50; T22 = 60;
Z=c1*exp(-TE/T21) + c2*exp(-TE/T22);
ZMAX=max(max(abs(Z)));
STDnoise=ZMAX*(1/SNR);

TT21_1D = zeros(nruns,1);
TT22_1D = zeros(nruns,1);
c1_1D = zeros(nruns,1);
c2_1D = zeros(nruns,1);

surfit = @(B,TE) B(3)*exp(-TE/B(1))+ B(4)*exp(-TE/B(2));
parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
% noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T22+10*rand() rand() rand()], TE, ZN, [],  [],opt);
TT21_1D(kk)=BN(1);TT22_1D(kk)=BN(2);c1_1D(kk) = BN(3); c2_1D(kk) = BN(4);
end
for i = 1:nruns
    if TT21_1D(i) > TT22_1D(i)
        a = TT21_1D(i);
        b = TT22_1D(i);
        TT21_1D(i) = b;
        TT22_1D(i) = a;
        
        c = c1_1D(i);
        d = c2_1D(i);
        c1_1D(i) = d;
        c2_1D(i) = c;
    end
end
sd_c1(1,1) = std(c1_1D);
sd_c2(1,1) = std(c2_1D);
sd_T21(1,1) = std(TT21_1D);
sd_T22(1,1) = std(TT22_1D);

%% plot
% fprintf(['SD(c1) = ', num2str(std(c1_1D)), ' SD(c2) = ', num2str(std(c1_1D)),...
%     ' SD(T21) = ', num2str(std(TT21_1D)), ' SD(T21) = '])

h1= figure;
tiledlayout(1,2)

ax1 = nexttile;
histogram(TT21_1D, nbins_T2)
hold on
histogram(TT22_1D, nbins_T2)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_1D),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_1D),2)]},'FontSize',30, 'Location', 'NorthEast') 
% title( {['T_{2,1}=' num2str(T21), ', T_{2,2}=' num2str(T22)]},'FontSize',30); 

ax2 = nexttile;
histogram(c1_1D, nbins_c)
hold on
histogram(c2_1D, nbins_c)
xlabel('c_1 or c_2','FontSize',30,'FontWeight','bold') 
legend({['c_1', ', SD = ', num2str(std(c1_1D),2)],['c_2', ', SD = ', num2str(std(c2_1D),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['1D: T_{2,1} = ' num2str(T21), ', T_{2,2} = ' num2str(T22)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1881         313        1284         401])

%% T21 = 60, T22 = 30

T21 = 47.5; T22 = 60;
Z=c1*exp(-TE/T21) + c2*exp(-TE/T22);
ZMAX=max(max(abs(Z)));
STDnoise=ZMAX*(1/SNR);

TT21_1D = zeros(nruns,1);
TT22_1D = zeros(nruns,1);
c1_1D = zeros(nruns,1);
c2_1D = zeros(nruns,1);

surfit = @(B,TE) B(3)*exp(-TE/B(1))+ B(4)*exp(-TE/B(2));
parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
% noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T22+10*rand() rand() rand()], TE, ZN, [],  [],opt);
TT21_1D(kk)=BN(1);TT22_1D(kk)=BN(2);c1_1D(kk) = BN(3); c2_1D(kk) = BN(4);
end
for i = 1:nruns
    if TT21_1D(i) > TT22_1D(i)
        a = TT21_1D(i);
        b = TT22_1D(i);
        TT21_1D(i) = b;
        TT22_1D(i) = a;
        
        c = c1_1D(i);
        d = c2_1D(i);
        c1_1D(i) = d;
        c2_1D(i) = c;
    end
end

%% plot
h2 = figure;
tiledlayout(1,2)

ax3 = nexttile;
histogram(TT21_1D, nbins_T2)
hold on
histogram(TT22_1D, nbins_T2)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_1D),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_1D),2)]},'FontSize',30, 'Location', 'NorthEast') 
% title( {['T_{2,1}=' num2str(T21), ', T_{2,2}=' num2str(T22)]},'FontSize',30); 

ax4 = nexttile;
histogram(c1_1D, nbins_c)
hold on
histogram(c2_1D, nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold') 
lgnd = legend({['c_1', ', SD = ', num2str(std(c1_1D),2)],['c_2', ', SD = ', num2str(std(c2_1D),2)]},'FontSize',30, 'Location', 'North') 
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
sgtitle( {['1D: T_{2,1} = ' num2str(T21), ', T_{2,2} = ' num2str(T22)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1881         313        1284         401])

%% T21 = 60, T22 = 30

T21 = 45; T22 = 60;
Z=c1*exp(-TE/T21) + c2*exp(-TE/T22);
ZMAX=max(max(abs(Z)));
STDnoise=ZMAX*(1/SNR);

TT21_1D = zeros(nruns,1);
TT22_1D = zeros(nruns,1);
c1_1D = zeros(nruns,1);
c2_1D = zeros(nruns,1);

surfit = @(B,TE) B(3)*exp(-TE/B(1))+ B(4)*exp(-TE/B(2));
parfor kk=1:nruns
noise = normrnd(0,STDnoise,size(Z));
% noise = normrnd(0,STDnoise,size(Z));
ZN=Z+noise;
BN = lsqcurvefit(surfit, [T21+10*rand() T22+10*rand() rand() rand()], TE, ZN, [],  [],opt);
TT21_1D(kk)=BN(1);TT22_1D(kk)=BN(2);c1_1D(kk) = BN(3); c2_1D(kk) = BN(4);
end
for i = 1:nruns
    if TT21_1D(i) > TT22_1D(i)
        a = TT21_1D(i);
        b = TT22_1D(i);
        TT21_1D(i) = b;
        TT22_1D(i) = a;
        
        c = c1_1D(i);
        d = c2_1D(i);
        c1_1D(i) = d;
        c2_1D(i) = c;
    end
end

%% plot
h3 = figure;
tiledlayout(1,2)

ax5 = nexttile;
histogram(TT21_1D, nbins_T2)
hold on
histogram(TT22_1D, nbins_T2)
xlabel('T_{2,1} or T_{2,2}','FontSize',40,'FontWeight','bold')
legend({['T_{2,1}', ', SD = ', num2str(std(TT21_1D),2)],['T_{2,2}', ', SD = ', num2str(std(TT22_1D),2)]},'FontSize',30, 'Location', 'NorthEast') 
% title( {['T_{2,1}=' num2str(T21), ', T_{2,2}=' num2str(T22)]},'FontSize',30); 

ax6 = nexttile;
histogram(c1_1D, nbins_c)
hold on
histogram(c2_1D, nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold') 
lgnd = legend({['c_1', ', SD = ', num2str(std(c1_1D),2)],['c_2', ', SD = ', num2str(std(c2_1D),2)]},'FontSize',30, 'Location', 'North') 
lgnd.BoxFace.ColorType='truecoloralpha';
lgnd.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
sgtitle( {['1D: T_{2,1} = ' num2str(T21), ', T_{2,2} = ' num2str(T22)]},'FontSize',40,'FontWeight','bold'); 
set(gcf,'position',[1881         313        1284         401])

%%
linkaxes([ax1,ax3,ax5],'xy')
linkaxes([ax2,ax4,ax6],'xy')

%%
saveas(h1,'MC_1D_slices_A','epsc')
saveas(h2,'MC_1D_slices_B','epsc')
saveas(h3,'MC_1D_slices_C','epsc')
