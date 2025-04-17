function par = Jacobian(x,u,par)

% This function computes a second order approximation of the Jacobian.

A = speye(par.n);
x_k = zeros(par.n,par.m);
    for k = 1:par.K-1
         A_bar = par.A;

    for r = 1:par.m
         A_bar = A_bar + u(r,par.K-k)*par.B(:,(r-1)*par.n+1:r*par.n);
         x_k(:,r) = par.B(:,(r-1)*par.n+1:r*par.n)*x(1:par.n,par.K-k);
    end

    A_bold = expm(par.dt*A_bar);
    B_bold = par.dt*A_bold*x_k;

    par.J(:,par.m*(par.K-k-1)+1:par.m*(par.K-k)) = A*B_bold;
    A = A*A_bold;
    end
end