%RGS 2/25/19
%wide/narrow
clear;clc;
close all;
%nruns is the number of fits for each value of T1,2.  nfits=500 or so is a
%good value for plots that show the effects
nruns=1000;
c1=0.3;c2=1-c1;
% T21=70;T22=85;
T21=45;T22=60;
TEmin=8;TEdel=8;TEmax=512;
%next line defines the range of the input simulation variable T12. Fits are
%done nruns times (each time with a different noise reaalization) for each value of T12
%A CHANGE
T12min=200; T12max=1800;delT12=20;

TI = 50:200:5000;
TE=TEmin:TEdel:TEmax;
%A CHANGE
T11=1000;
opt = optimoptions('lsqcurvefit','Display','off','Algorithm','levenberg-marquardt');

[X,Y] = meshgrid(TE,TI);

XY(:,:,1) = X;
XY(:,:,2) = Y;

surfit = @(B,XY) B(5)*exp(-TE/B(1)) .* (1-2.*exp((-TI/B(2))))'+ B(6)*exp(-TE/B(3)) .* (1-2.*exp((-TI/B(4))))';
BN=zeros(1,6);

SNR=400;


%calculating SNR for 1st values of parameters, and holding noise level
%constant for all subsequent runs, to correspond to experiments.
Z=c1*exp(-TE/T21).*(1-2.*exp(-TI/T11))' + c2*exp(-TE/T22).*(1-2.*exp(-TI/T12min))';
ZMAX=max(max(abs(Z)));  
STDnoise=ZMAX*(1/SNR);
% STDnoise=0.05
%% run MC inversion

T12_vec = T12min:delT12:T12max;
N = length(T12_vec);
TT21 = zeros(N,nruns);
TT11 = zeros(N,nruns);
TT22 = zeros(N,nruns);
TT12 = zeros(N,nruns);
cc1 = zeros(N,nruns);
cc2 = zeros(N,nruns);
for n = 1:N
    T12 = T12_vec(n);
    Z=c1*exp(-TE/T21).*(1-2.*exp(-TI/T11))' + c2*exp(-TE/T22).*(1-2.*exp(-TI/T12))';
        parfor kk=1:nruns
        %the next line gives the size of the noise
        noise = normrnd(0,STDnoise,size(Z));
        %the next line gives the noisy signal for the particular value of
        %T12 in this loop; remember there are "nruns" different times this
        %is done, each with its own fit
        ZN=Z+noise;
        %The next line is the fit of the noisy data for the 6 parameters,
        %T21, T11,T22,T12,c1, and c2
        % Initial guess to be the true solution
%         BN = lsqcurvefit(surfit, [50*rand() 50*rand() 500*rand() 500*rand() rand() rand()], XY, ZN,[],[],opt);
        BN = lsqcurvefit(surfit, [T21+10*rand() T11+20*rand() T22+10*rand() T12+20*rand() c1+rand() c2+rand()], XY, ZN,[],[],opt);
        
        if BN(1)>BN(3)
            a = BN(5);
            BN(5) = BN(6);
            BN(6) = a;
            
            b = BN(1);
            BN(1) = BN(3);
            BN(3) = b;
        end
         
        
        TT21(n,kk)=BN(1);TT11(n,kk)=BN(2);TT22(n,kk)=BN(3);
        TT12(n,kk)=BN(4);cc1(n,kk)=BN(5);cc2(n,kk)=BN(6);
        end
end

%% plot
T12_vec_plot1= T12_vec' * ones(1,nruns);
T12_vec_plot = T12_vec_plot1(:);

% T21
h1 = figure;
subplot(2,1,1) 
scatter(T12_vec_plot,TT21(:),'*')
ylabel('T_{2,1}^*','FontSize',36,'FontWeight','bold')
% xlabel('Underlying T_{1,2}','FontSize',36,'FontWeight','bold')

var_TT21 = zeros(size(TT21,1),1);
for i = 1:length(var_TT21)
    var_TT21(i) = var(TT21(i,:));
end
subplot(2,1,2)
plot(T12min:delT12:T12max, var_TT21,'LineWidth',3);
ylabel('SD(T_{2,1}^*)','FontSize',36,'FontWeight','bold')
xlabel('Underlying T_{1,2}','FontSize',36,'FontWeight','bold')

saveas(h1,'MC_2D_T2A_std','epsc')
%%
% T22
h2 = figure;
subplot(2,1,1) 
scatter(T12_vec_plot,TT22(:),'*')
ylabel('T_{2,2}^*','FontSize',36,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',36,'FontWeight','bold')

var_TT22 = zeros(size(TT22,1),1);
for i = 1:length(var_TT22)
    var_TT22(i) = var(TT22(i,:));
end
subplot(2,1,2)
plot(T12min:delT12:T12max, var_TT22,'LineWidth',3);
ylabel('SD(T_{2,2}^*)','FontSize',36,'FontWeight','bold')
xlabel('Underlying T_{1,2}','FontSize',36,'FontWeight','bold')
saveas(h2,'MC_2D_T2B_std','epsc')

%%
% % T12
% figure
% subplot(2,1,1) 
% 
% scatter(T12_vec_plot,TT12(:), '*')
% 
% ylabel('T_{1,2} Fit Values ','FontSize',26,'FontWeight','bold')
% xlabel('  T_{1,B} Values ','FontSize',26,'FontWeight','bold')
% 
% var_TT12 = zeros(size(TT12,1),1);
% for i = 1:length(var_TT12)
%     var_TT12(i) = var(TT12(i,:));
% end
% subplot(2,1,2)
% plot(T12min:delT12:T12max, var_TT12,'LineWidth',3);
% ylabel('SD(T_{1,2})','FontSize',26,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',26,'FontWeight','bold')
% 
% %%
% 
% % T11
% figure
% subplot(2,1,1) 
% 
% scatter(T12_vec_plot,TT11(:), '*')
% 
% ylabel('T_{1,1} Fit Values ','FontSize',26,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',26,'FontWeight','bold')
% var_TT11 = zeros(size(TT11,1),1);
% for i = 1:length(var_TT11)
%     var_TT11(i) = var(TT11(i,:));
% end
% subplot(2,1,2)
% plot(T12min:delT12:T12max, var_TT11,'LineWidth',3);
% ylabel('SD(T_{1,1})','FontSize',26,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',26,'FontWeight','bold')

%%
% c1
h3 = figure;
subplot(2,1,1) 

scatter(T12_vec_plot,cc1(:), '*')
ylabel('c_1^*','FontSize',36,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',36,'FontWeight','bold')
var_cc1 = zeros(size(cc1,1),1);
for i = 1:length(var_cc1)
    var_cc1(i) = var(cc1(i,:));
end
subplot(2,1,2)
plot(T12min:delT12:T12max, var_cc1,'LineWidth',3);
ylabel('SD(c_1^*)','FontSize',36,'FontWeight','bold')
xlabel('Underlying T_{1,2}','FontSize',36,'FontWeight','bold')
saveas(h3,'MC_2D_cA_std','epsc')

%%
% c2
h4 = figure;
subplot(2,1,1) 

scatter(T12_vec_plot,cc2(:), '*')
ylabel('c_2^*','FontSize',36,'FontWeight','bold')
% xlabel('T_{1,2}','FontSize',36,'FontWeight','bold')
var_cc2 = zeros(size(cc2,1),1);
for i = 1:length(var_cc2)
    var_cc2(i) = var(cc2(i,:));
end
subplot(2,1,2)
plot(T12min:delT12:T12max, var_cc2,'LineWidth',3);
ylabel('SD(c_2^*)','FontSize',36,'FontWeight','bold')
xlabel('Underlying T_{1,2}','FontSize',36,'FontWeight','bold')
saveas(h4,'MC_2D_cB_std','epsc')
