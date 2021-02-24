function D = Derivatives_3D(TE,tau,b,P_3D)
%10/10/2018
% these are the derivatives of the model; there is no specific "Jacobian"
% since there are three "time" variables, so this is "by hand"
% note:  P_3D=[x1 x2 x3 x4 x5 x6 x7 x8]
%x1 is c1; x2 is T2,1; x3 is T1,1; x4 is ADC,1
%x5 is c2; x6 is T2,2; x7 is T1,2; x8 is ADC,2
% model, written in the format for T2,T1,D experiments 
% with vectors of variables 
% t_i (this is directly detected)
% tau_j (this is indirectly detected, but looped over to see the errors as a function of T1)
% b_k (this is indirectly detected, but looped over to see the errors as a function of ADC) 
% 2.5 * 10^-3 mm^2/sec water at room temperature; tissue range ~1.0 - 2.5

% s(x,t_i,tau_j,b_k)=x1*exp(-t_i/x2)*(1-2exp(-tau_j/x3))*exp(-b_k*x4)+x5*exp(-t_i/x6)*(1-2exp(-tau_j/x7))*exp(-b_k*x8)

%note:  no offset in the above model; this can be added if desired

nparms=length(P_3D);
nTE=length(TE);
ntau=length(tau);
nb=length(b);

D=zeros(nTE*ntau*nb, nparms);
%rename the direct dimension variable from TE to t; the 1st indirect dimension
%time variable is tau, the 2nd indirect dimension variable is b
t=TE;

x1=P_3D(1);
x2=P_3D(2);
x3=P_3D(3);
x4=P_3D(4);
x5=P_3D(5);
x6=P_3D(6);
x7=P_3D(7);
x8=P_3D(8);

% note:  P_3D=[x1 x2 x3 x4 x5 x6 x7 x8]
%x1 is c1; x2 is T2,1; x3 is T1,1; x4 is ADC,1
%x5 is c2; x6 is T2,2; x7 is T1,2; x8 is ADC,2

ix=0;
for i=1:nTE
    for j=1:ntau
        for k=1:nb
            
        ix=ix+1;
       
        D(ix,1)=exp(-t(i)/x2) * (1-2*exp(-tau(j)/x3)) * exp(-b(k)*x4);                   %Derivative wrt x1 
        D(ix,2)=x1*t(i)/(x2^2)*exp(-t(i)/x2)*(1-2*exp(-tau(j)/x3)) * exp(-b(k)*x4);      %Derivative wrt x2 
        D(ix,3)= -2*x1*(tau(j)/(x3^2))*exp(-t(i)/x2)*(1-2*exp(-tau(j)/x3))* exp(-b(k)*x4); %Derivative wrt x3 
        D(ix,4)=x1*exp(-t(i)/x2)*(1-2*exp(-tau(j)/x3)) * -b(k)*exp(-b(k)*x4);         %Derivative wrt x4 
        
        D(ix,5)=exp(-t(i)/x6) * (1-2*exp(-tau(j)/x7))*exp(-b(k)*x8);                   %Derivative wrt x5
        D(ix,6)=x5*t(i)/(x6^2)*exp(-t(i)/x6)*(1-2*exp(-tau(j)/x7))*exp(-b(k)*x8);    %Derivative wrt x6 
        D(ix,7)= -2*x5*(tau(j)/(x7^2))*exp(-t(i)/x6)*(1-2*exp(-tau(j)/x7))*exp(-b(k)*x8); %Derivative wrt x7 
        D(ix,8)=x5*exp(-t(i)/x6)*(1-2*exp(-tau(j)/x7))* -b(k)*exp(-b(k)*x8);        %Derivative wrt x8
        
        end
    end
end
end
