function D = T2Derivatives_biexp_no_const(TE,P_biexp)
%these are the derivatives of the model; the (-) sign to form the Jacobian
%will be added in the calling code, although this doesn't make any
%numerical difference since the cov matrix is the product of two Jacobian
%terms
%note:  P_biexp=[c2 c3 tau1 tau2]

%either group of statements below works and gives the same answers.  Evidently, MATLAB understand that
%the first group of statements produces column vectors, as desired, because
%the LHS is defined as a col vector, even though the RHS expressions are
%row vectors
%%model:  M(x,TE)=c1*exp(-TE/tau1)+c2*exp(-TE/tau2)
%note:  P_biexp=[c1 c2 tau1 tau2]
D(:,1)=exp(-TE/P_biexp(3));        %Derivative wrt c1
D(:,2)=exp(-TE/P_biexp(4));        %Derivative wrt c2
D(:,3)=P_biexp(1)*TE.*exp(-TE/P_biexp(3))./(P_biexp(3)^2);    %Derivative wrt tau1
D(:,4)=P_biexp(2)*TE.*exp(-TE/P_biexp(4))./(P_biexp(4)^2);    %Derivative wrt tau2


end
