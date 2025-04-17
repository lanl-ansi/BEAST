function [error,x] = simulate_base(par,u,Schrodinger,a,b)

% This function simulates the state of the original system corresponding to
% the control action computed with state transfer control algorithm 
% (or minimal energy if min_energy = true). If par.Na ~= 1 or par.Nb ~= 1,
% then the simulation is repeated for different values of the variable
% parameters alpha and beta defiend by a and b.

error = zeros(length(a),length(b));
x = [];

for aa = 1:length(a)
    for bb = 1:length(b)
        y(:,1) = par.X0(1:par.n_base);
    
        for k=1:par.K-1    
            A = a(aa)*par.A_base;                           
            for r = 1:par.m
                A = A + b(bb)*u(r,k)*par.B_base(:,(r-1)*par.n_base+1:r*par.n_base);
            end

            y(1:par.n_base,k+1) = expm(par.dt*A)*y(1:par.n_base,k);   
        end

        if Schrodinger == true
            error(aa,bb) = norm((y(1:par.n_base/2,end).^2 + y(par.n_base/2+1:end,end).^2).^(1/2)-full(par.XT(1:par.n_base/2)))^2;
            x = [x; y(1:par.n_base/2,:).^2 + y(par.n_base/2+1:end,:).^2];
        else
            error(aa,bb) = norm(y(:,end)-par.XT(1:par.n_base))^2;
            x = [x; y];
        end
    end
end

end