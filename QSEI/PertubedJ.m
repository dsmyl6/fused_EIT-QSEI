function [J] = PertubedJ_complex(nel,theta,nu,th,constraint,Tri,g,cond,force)
%%% Function to obtain Jacobian using complex number pertubation
[z,~]=FMDL(theta,nu,th,constraint,Tri,g,cond,force); % function value
n=numel(theta);                     % size of independent
m=numel(z);                     % size of dependent
J=zeros(m,n);                   % allocate memory for the Jacobian matrix
h=n*eps;                        % differentiation step size
for k=1:n                       % loop for each independent variable 
    x1=theta;                       % reference point
    x1(k)=x1(k)+h*i;            % increment in kth independent variable
    J(:,k)=imag(FMDL(x1,nu,th,constraint,Tri,g,cond,force))/h;     % complex step differentiation 
end