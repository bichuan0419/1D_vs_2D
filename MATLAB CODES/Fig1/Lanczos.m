%% Lanczos example
%% A variation of Lanczos example: two exponential case
clc;clear;close all
% Lanczos example
x = linspace(0,400,24)';
ob1 = 0.2*exp(-x/20) + 0.8*exp(-x/150);

% find another two-component expoential decays 
fun = @(y) (y(1)*exp(-x/40) + y(2)*exp(-x/y(3)) - ob1)'*(y(1)*exp(-x/40) + y(2)*exp(-x/y(3)) - ob1);
fun_1 = @(y) y(1)*exp(-x/40) + y(2)*exp(-x/y(3)) ;
y0 = [0.4,0.6,60];
y_approx = fmincon(fun,y0,[],[]);
ob3 = fun_1(y_approx);

ax1 = figure;
hold on
plot(x,ob1,'-.or','LineWidth',3);
% plot(x,ob2,'-.*b','LineWidth',1.5);
plot(x,ob3,'-.xg','LineWidth',3);
xlabel('t','FontSize',30,'FontWeight','bold')
ylabel('Amplitude','FontSize',30,'FontWeight','bold')
legend({'Model 1', 'Model 2'},'FontSize',24)

saveas(ax1,'Lanczos','epsc')
