function x = simulate(par,u)

% This function simulates the state of the Legendre expansion coefficients
% corresponding to the supplied control action for each iteration of the
% sequential quadratic program.

x(:,1) = par.X0;

for k = 1:par.K-1    
    A = par.A;                           
      
    for r = 1:par.m
          A = A + u(r,k)*par.B(:,(r-1)*par.n+1:r*par.n);
    end
     
    x(1:par.n,k+1) = expm(par.dt*A)*x(1:par.n,k);                    
end

end
