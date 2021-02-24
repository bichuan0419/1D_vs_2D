clear;clc;
close all;
nruns = 1000;
SNR = 1;
nbins_nu = 20;
nbins_c = 20;
nbins_mu = 20;

%% 1D
L = 20^2;             % Length of signal
L_A = 10^2;           % Length of the second dimension
t = linspace(0,1,L);        % Time vector
tau = linspace(0,1,L_A);        % Time vector
T = t(2) - t(1);
Fs = 1/T;
Fs_A = tau(2) - tau(1);
c1 = 0.7;
c2 = 0.3;
nu1 = 30;
nu2 = 40;
%% 1D
% define forward problem function
% S_1D = c1*sin(2*pi*nu1*t) + c2*sin(2*pi*nu2*t);
S_1D = c1*exp(2*pi*1i*nu1*t) + c2*exp(2*pi*1i*nu2*t);
std_1D = max(abs(S_1D))/SNR;
c1_vec_1D = zeros(nruns,1);
c2_vec_1D = zeros(nruns,1);
nu1_vec_1D = zeros(nruns,1);
nu2_vec_1D = zeros(nruns,1);
f = Fs*(0:(L/2))/L;
%%
for i =1 : nruns
    
    SN = S_1D + normrnd(0,std_1D/sqrt(L_A),size(S_1D));
    Y = fft(SN);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    [max_vec, id_vec] = maxk(P1,2);
    max_vec = max_vec/sum(max_vec);
    c1_vec_1D(i) = max_vec(1);
    c2_vec_1D(i) = max_vec(2);
    nu1_vec_1D(i) = f(id_vec(1));
    nu2_vec_1D(i) = f(id_vec(2));
end

h1 = figure;
tiledlayout(1,2)

ax1 = nexttile;
histogram(nu1_vec_1D,nbins_nu)
hold on
histogram(nu2_vec_1D,nbins_nu)
xlabel('\nu_1 or \nu_2','FontSize',40,'FontWeight','bold')
legend({['\nu_1', ', \sigma = ', num2str(std(nu1_vec_1D),2)],['\nu_2', ', \sigma = ', num2str(std(nu2_vec_1D),2)]},'FontSize',30, 'Location', 'North') 
% title({['1D: \nu_1 = ', num2str(nu1), ', \nu_2 = ', num2str(nu2)]},'FontSize',22)


ax2 = nexttile;
histogram(c1_vec_1D,nbins_c)
hold on
histogram(c2_vec_1D,nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold')
legend({['c_1', ', \sigma = ', num2str(std(c1_vec_1D),2)],['c_2', ', \sigma = ', num2str(std(c2_vec_1D),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['1D: \nu_1 = ' num2str(nu1), ', \nu_2 = ' num2str(nu2)]},'FontSize',40,'FontWeight','bold'); 

set(gcf,'position',[1845         159        1284         348])

%% 2D
mu1 = 200;
mu2 = 200;
S_2D = zeros(L,L_A);
for i = 1:L
    for j = 1:L_A
%         S_2D(i,j) = c1*sin(2*pi*(nu1*t(i) + mu1*tau(j))) + c2*sin(2*pi*(nu2*t(i) + mu2*tau(j)));

         S_2D(i,j) = c1*exp(2*pi*1i*(nu1*t(i) + mu1*tau(j))) + c2*exp(2*pi*1i*(nu2*t(i) + mu2*tau(j)));
    end
end
std_2D = max(max(abs(S_2D)))/SNR;
c1_vec_2D = zeros(nruns,1);
c2_vec_2D = zeros(nruns,1);
nu1_vec_2D = zeros(nruns,1);
nu2_vec_2D = zeros(nruns,1);
mu1_vec_2D = zeros(nruns,1);
mu2_vec_2D = zeros(nruns,1);

f = Fs*(0:(L/2))/L;
f_A = Fs_A*(0:(L_A/2))/L_A;
%%
for i =1 : nruns
    
    SN = S_2D + normrnd(0,std_2D,size(S_2D));
    Y = fft2(SN);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1,1:L_A/2+1);
    temp = P1;
    P1(:,2:end-1) = 2*temp(:,2:end-1);
    P1(2:end-1,:) = 2*temp(2:end-1,:);

    [max_vec, id_rec] = maxk(P1(:),2);
    c_rec = max_vec/sum(max_vec);
    
    [id_nu1_vec_2D , id_mu1_vec_2D] = ind2sub(size(P1),id_rec(1));
    nu1_vec_2D(i) = f(id_nu1_vec_2D);
    mu1_vec_2D(i) = f_A(id_mu1_vec_2D);
    
    [id_nu2_vec_2D , id_mu2_vec_2D] = ind2sub(size(P1),id_rec(2));
    nu2_vec_2D(i) = f(id_nu2_vec_2D);
    mu2_vec_2D(i) = f_A(id_mu2_vec_2D);
    
    c1_vec_2D(i) = c_rec(1);
    c2_vec_2D(i) = c_rec(2);
end
%%
h2 = figure;
tiledlayout(1,2)

ax3 = nexttile;
histogram(nu1_vec_2D,nbins_nu)
hold on
histogram(nu2_vec_2D,nbins_nu)
xlabel('\nu_1 or \nu_2','FontSize',40,'FontWeight','bold')
legend({['\nu_1', ', \sigma = ', num2str(std(nu1_vec_2D),2)],['\nu_2', ', \sigma = ', num2str(std(nu2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
% title({['1D: \nu_1 = ', num2str(nu1), ', \nu_2 = ', num2str(nu2)]},'FontSize',22)


ax4 = nexttile;
histogram(c1_vec_2D,nbins_c)
hold on
histogram(c2_vec_2D,nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold')
legend({['c_1', ', \sigma = ', num2str(std(c1_vec_2D),2)],['c_2', ', \sigma = ', num2str(std(c2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: \mu_1 = ' num2str(mu1), ', \mu_2 = ' num2str(mu2)]},'FontSize',40,'FontWeight','bold'); 

set(gcf,'position',[1845         159        1284         348])
%% 2D
mu1 = 200;
mu2 = 100;
S_2D = zeros(L,L_A);
for i = 1:L
    for j = 1:L_A
        S_2D(i,j) = c1*sin(2*pi*(nu1*t(i) + mu1*tau(j))) + c2*sin(2*pi*(nu2*t(i) + mu2*tau(j)));
    end
end
std_2D = max(max(S_2D))/SNR;
c1_vec_2D = zeros(nruns,1);
c2_vec_2D = zeros(nruns,1);
nu1_vec_2D = zeros(nruns,1);
nu2_vec_2D = zeros(nruns,1);
mu1_vec_2D = zeros(nruns,1);
mu2_vec_2D = zeros(nruns,1);

f = Fs*(0:(L/2))/L;
f_A = Fs_A*(0:(L_A/2))/L_A;
%%
for i =1 : nruns
    
    SN = S_2D + normrnd(0,std_2D,size(S_2D));
    Y = fft2(SN);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1,1:L_A/2+1);
    temp = P1;
    P1(:,2:end-1) = 2*temp(:,2:end-1);
    P1(2:end-1,:) = 2*temp(2:end-1,:);

    [max_vec, id_rec] = maxk(P1(:),2);
    c_rec = max_vec/sum(max_vec);
    
    [id_nu1_vec_2D , id_mu1_vec_2D] = ind2sub(size(P1),id_rec(1));
    nu1_vec_2D(i) = f(id_nu1_vec_2D);
    mu1_vec_2D(i) = f_A(id_mu1_vec_2D);
    
    [id_nu2_vec_2D , id_mu2_vec_2D] = ind2sub(size(P1),id_rec(2));
    nu2_vec_2D(i) = f(id_nu2_vec_2D);
    mu2_vec_2D(i) = f_A(id_mu2_vec_2D);
    
    c1_vec_2D(i) = c_rec(1);
    c2_vec_2D(i) = c_rec(2);
end
%%
h3 = figure;
tiledlayout(1,2)

ax5 = nexttile;
histogram(nu1_vec_2D,nbins_nu)
hold on
histogram(nu2_vec_2D,nbins_nu)
xlabel('\nu_1 or \nu_2','FontSize',40,'FontWeight','bold')
legend({['\nu_1', ', \sigma = ', num2str(std(nu1_vec_2D),2)],['\nu_2', ', \sigma = ', num2str(std(nu2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
% title({['1D: \nu_1 = ', num2str(nu1), ', \nu_2 = ', num2str(nu2)]},'FontSize',22)


ax6 = nexttile;
histogram(c1_vec_2D,nbins_c)
hold on
histogram(c2_vec_2D,nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold')
legend({['c_1', ', \sigma = ', num2str(std(c1_vec_2D),2)],['c_2', ', \sigma = ', num2str(std(c2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: \mu_1 = ' num2str(mu1), ', \mu_2 = ' num2str(mu2)]},'FontSize',40,'FontWeight','bold'); 

set(gcf,'position',[1845         159        1284         348])
%% 2D
mu1 = 200;
mu2 = 50;
S_2D = zeros(L,L_A);
for i = 1:L
    for j = 1:L_A
        S_2D(i,j) = c1*sin(2*pi*(nu1*t(i) + mu1*tau(j))) + c2*sin(2*pi*(nu2*t(i) + mu2*tau(j)));
    end
end
std_2D = max(max(S_2D))/SNR;
c1_vec_2D = zeros(nruns,1);
c2_vec_2D = zeros(nruns,1);
nu1_vec_2D = zeros(nruns,1);
nu2_vec_2D = zeros(nruns,1);
mu1_vec_2D = zeros(nruns,1);
mu2_vec_2D = zeros(nruns,1);

f = Fs*(0:(L/2))/L;
f_A = Fs_A*(0:(L_A/2))/L_A;
%%
for i =1 : nruns
    
    SN = S_2D + normrnd(0,std_2D,size(S_2D));
    Y = fft2(SN);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1,1:L_A/2+1);
    temp = P1;
    P1(:,2:end-1) = 2*temp(:,2:end-1);
    P1(2:end-1,:) = 2*temp(2:end-1,:);

    [max_vec, id_rec] = maxk(P1(:),2);
    c_rec = max_vec/sum(max_vec);
    
    [id_nu1_vec_2D , id_mu1_vec_2D] = ind2sub(size(P1),id_rec(1));
    nu1_vec_2D(i) = f(id_nu1_vec_2D);
    mu1_vec_2D(i) = f_A(id_mu1_vec_2D);
    
    [id_nu2_vec_2D , id_mu2_vec_2D] = ind2sub(size(P1),id_rec(2));
    nu2_vec_2D(i) = f(id_nu2_vec_2D);
    mu2_vec_2D(i) = f_A(id_mu2_vec_2D);
    
    c1_vec_2D(i) = c_rec(1);
    c2_vec_2D(i) = c_rec(2);
end
%%
h4 = figure;
tiledlayout(1,2)

ax7 = nexttile;
histogram(nu1_vec_2D,nbins_nu)
hold on
histogram(nu2_vec_2D,nbins_nu)
xlabel('\nu_1 or \nu_2','FontSize',40,'FontWeight','bold')
legend({['\nu_1', ', \sigma = ', num2str(std(nu1_vec_2D),2)],['\nu_2', ', \sigma = ', num2str(std(nu2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
% title({['1D: \nu_1 = ', num2str(nu1), ', \nu_2 = ', num2str(nu2)]},'FontSize',22)


ax8 = nexttile;
histogram(c1_vec_2D,nbins_c)
hold on
histogram(c2_vec_2D,nbins_c)
xlabel('c_1 or c_2','FontSize',40,'FontWeight','bold')
legend({['c_1', ', \sigma = ', num2str(std(c1_vec_2D),2)],['c_2', ', \sigma = ', num2str(std(c2_vec_2D),2)]},'FontSize',30, 'Location', 'North') 
sgtitle( {['2D: \mu_1 = ' num2str(mu1), ', \mu_2 = ' num2str(mu2)]},'FontSize',40,'FontWeight','bold'); 

set(gcf,'position',[1845         159        1284         348])
% %% save fig
saveas(h1,'IFT_1D_vs_2D_A','epsc')
saveas(h2,'IFT_1D_vs_2D_B','epsc')
saveas(h3,'IFT_1D_vs_2D_C','epsc')
saveas(h4,'IFT_1D_vs_2D_D','epsc')
%%
std(c1_vec_1D)
std(c1_vec_2D)
std(nu1_vec_1D)
std(nu1_vec_2D)