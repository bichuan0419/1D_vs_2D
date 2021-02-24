function D = partial_derivatives_2D(TE,tau,P_2D)
%these are the derivatives of the model; there is no specific "Jacobian"
%since there are two time variables, so this is "by hand", but the same
%idea
%note:  P_2D=[x1 x2 x3 x4 x5 x6] = [c1 T2A T1A c2 T2B T1B]
%model:  s(x,t_i, tau_j)=x1*exp(-t_i/x2)*(1-2exp(-tau_j/x3))+x4*exp(-t_i/x5)*(1-2exp(-tau_j/x6))
%note:  no offset in the above model; this can be added if desired
nparms=length(P_2D);
nTE=length(TE);
ntau=length(tau);
D=zeros(nTE*ntau, nparms);
%rename the direct dimension variable from TE to t; the indirect dimension
%time variable is tau
t=TE;
x1=P_2D(1);
x2=P_2D(2);
x3=P_2D(3);
x4=P_2D(4);
x5=P_2D(5);
x6=P_2D(6);
ix=0;
for j=1:ntau
    for i=1:nTE

        ix=ix+1;
        D(ix,1)= exp(-t(i)/x2) * (1-2*exp(-tau(j)/x3));                  %Derivative wrt x1 
        D(ix,2)=x1*t(i)/(x2^2)*exp(-t(i)/x2)*(1-2*exp(-tau(j)/x3));      %Derivative wrt x2 
        D(ix,3)=-2*x1*tau(j)/(x3^2)*exp(-t(i)/x2)*(1-2*exp(-tau(j)/x3)); %Derivative wrt x3 
        D(ix,4)=exp(-t(i)/x5) * (1-2*exp(-tau(j)/x6));                   %Derivative wrt x4 
        D(ix,5)=x4*t(i)/(x5^2)*exp(-t(i)/x5)*(1-2*exp(-tau(j)/x6));  %Derivative wrt x5 
        D(ix,6)=-2*x4*tau(j)/(x6^2)*exp(-t(i)/x5)*(1-2*exp(-tau(j)/x6)); %Derivative wrt x6 
    end
end
%the following commented lines are from the biexponential code
% D(:,1)=ones(length(TE),1);         %Derivative wrt c1 
% D(:,2)=exp(-TE/P_biexp(4));        %Derivative wrt c2
% D(:,3)=exp(-TE/P_biexp(5));        %Derivative wrt c3
% D(:,4)=P_biexp(2)*TE.*exp(-TE/P_biexp(4))./(P_biexp(4)^2);    %Derivative wrt tau1
% D(:,5)=P_biexp(3)*TE.*exp(-TE/P_biexp(5))./(P_biexp(5)^2);    %Derivative wrt tau2
end
