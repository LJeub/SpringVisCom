function [x_new, grad_new, step_new,L_new]=BFGS_update(x, grad,step,L,c1,c2,gradient,energy)
% compute the BFGS update step

% Version: 1.2
% Date: Tue 13 May 2014 17:03:19 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk


%This maintains an approximation for the Hessian
%
% g=L\(-grad(:));
% d=L'\g;
% d=reshape(d,size(grad));
% 
% step_new=wolfe_step(x,d,step,c1,c2,gradient,energy);
% delta=step_new*d;
% x_new=x+delta;
% grad_new=gradient(x_new);
% B=L*L';
% gamma=grad_new(:)-grad(:);
% beta=sqrt(gamma'*delta(:)/(delta(:)'*B*delta(:)));
% J=L'+L'*(delta(:)*(gamma-beta*B*delta(:))')/(beta*delta(:)'*B*delta(:));
% [~,L_new]=qr(J);
% L_new=L_new';

%This maintains an approximation for the inverse Hessian

d=-L*grad(:);
d=reshape(d,size(grad));
step_new=wolfe_step(x,d,step,c1,c2,gradient,energy);
x_new=x+step_new*d;
grad_new=gradient(x_new);
s=x_new(:)-x(:);
y=grad_new(:)-grad(:);
r=1/(y'*s);
I=eye(size(L));
L_new=(I-r*s*y')*L*(I-r*y*s')+r*(s*s');

end
