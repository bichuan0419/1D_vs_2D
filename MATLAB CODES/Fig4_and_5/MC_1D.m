%% 1D Monte-Carlo
clc;clear;close all
%% Parameter settings
% c1, c2=populations; tau1, tau2=time constants;
% TE=echo times, which are the sampling times
% model:  M(x,TE)=c1*exp(-TE/tau1)+c2*exp(-TE/tau2)
c1=0.7;
c2=(1-c1);
T21=60;

tau2min=30;
tau2max=100;
% del=(tau2max-tau2min)/100.;
del=1;
T22_vec=(tau2min:del:tau2max)';
nruns=length(T22_vec);
n_rlzn = 1000;
TE=(8:8:512)';


%% Create Objective Function; this is the fitting model: 

opt = optimoptions(@lsqcurvefit,'Display','off','Algorithm','Levenberg-Marquardt');
surfit = @(B,TE) B(1)*exp(-TE/B(3))+ B(2)*exp(-TE/B(4));
% surfit = @(B,TE) c1*exp(-TE/B(1))+ c2*exp(-TE/B(2));
%%
%defining the number of sub-plots

plotrows=2;
plotcols=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1_mtx=zeros(nruns,n_rlzn);
c2_mtx=zeros(nruns,n_rlzn);
T21_mtx=zeros(nruns,n_rlzn);
T22_mtx=zeros(nruns,n_rlzn);

parfor i = 1:nruns
    T22 = T22_vec(i);
    Z=c1*exp(-TE/T21) + c2*exp(-TE/T22);
    ZMAX=max(max(abs(Z)));
    SNR = 10000;
    STDnoise=ZMAX*(1/SNR);
    for j = 1:n_rlzn
        noise = normrnd(0,STDnoise,size(Z));
        ZN = Z + noise;
        BN = lsqcurvefit(surfit,[c1+rand() c2+rand() T21+10*rand() T22+10*rand()], TE, ZN,[],[],opt);
%         BN = lsqcurvefit(surfit,[abs(rand()) abs(rand()) abs(40*rand()) abs(100*rand())], TE, ZN,[],[],opt);

        if BN(3) > BN(4)
            a = BN(1);
            BN(1) = BN(2);
            BN(2) = a;
            
            b = BN(3);
            BN(3) = BN(4);
            BN(4) = b;
        end
%         if BN(1) >0 && BN(2) > 0 && BN(3)< 10000 && BN(4) < 10000
            c1_mtx(i,j) = BN(1); c2_mtx(i,j) = BN(2);
            T21_mtx(i,j) = BN(3); T22_mtx(i,j) = BN(4);
%         end
    end
end

%% plot wide-narrow scatter plots
% c1 variance
T22_vec_plot= T22_vec * ones(1,n_rlzn);
T22_vec_plot = T22_vec_plot(:);

%%
% swap ground truth
true_T22_vec = zeros(length(T22_vec),1);
for i = 1:length(T22_vec)
    if T22_vec(i) < T21
        true_T22_vec(i) = T21;
    else
        true_T22_vec(i) = T22_vec(i);
    end
end

% T22 variance
h4 = figure;
subplot(2,1,1) 
scatter(T22_vec_plot, abs(T22_mtx(:)),'*')
ylabel('T_{2,2}^*','FontSize',26,'FontWeight','bold')
set(gca,'yscale','log')
ylim([10 10^7])
yticks([10^1 10^2 10^3 10^4 10^5 10^6 10^7])
std_T22 = zeros(nruns,1);
for i = 1:nruns
    std_T22(i) = std(T22_mtx(i,:));
end
hold on
plot(T22_vec, true_T22_vec, 'LineWidth',2)
legend({'Fitted Value', 'Ground Truth'},'FontSize',20)
% ylim([0 500])

subplot(2,1,2)
semilogy(T22_vec, std_T22,'LineWidth',3);
xlabel('Underlying T_{2,2}','FontSize',26,'FontWeight','bold')
ylabel('SD(T_{2,2}^*)','FontSize',26,'FontWeight','bold')
yticks([10^0 10^1 10^2 10^3 10^4 10^5])
saveas(h4,'MC_1D_T22','epsc')

