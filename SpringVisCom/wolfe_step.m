function alpha=wolfe_step(x,dir,alpha_0,c_1,c_2,gradient,energy)
% find step size satisfying Wolfe conditions with parameters c_1, c_2

% Version: 1.2
% Date: Tue 13 May 2014 17:03:19 BST
% Author: Lucas G. S. Jeub
% Email: jeub@maths.ox.ac.uk

alpha=alpha_0;
alpha_low=0;
alpha_high=0;

e_0=energy(x);
grad_0=gradient(x);

ll=grad_0(:)'*dir(:);
cond1=energy(x+alpha*dir)>e_0+c_1*alpha*ll;
cond2=dphi(alpha)<c_2*ll;

while(cond1)||(cond2)
    update=false;
    if alpha-alpha_low<eps;
        error('stalled')
    end
    
    %alpha too large
    if cond1
        update=true;
    end
    
    while(cond1)
        alpha_high=alpha;
        alpha=(alpha_low+alpha_high)/2;
        cond1=energy(x+alpha*dir)>e_0+c_1*alpha*ll;
    end
    
    if update
        cond2=dphi(alpha)<c_2*ll;
    end
    %alpha too small
    if(cond2)
        alpha_low=alpha;
        if alpha_high==0
            alpha=2*alpha_low;
        else
            alpha=(alpha_low+alpha_high)/2;
        end
        cond1=energy(x+alpha*dir)>e_0+c_1*alpha*ll;
        cond2=dphi(alpha)<c_2*ll;
    end
end

    function d_phi=dphi(alpha)
        d_phi=gradient(x+alpha*dir);
        d_phi=d_phi(:)'*dir(:);
    end
end
