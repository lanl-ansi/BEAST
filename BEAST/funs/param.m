function par = param(par,Schrodinger)

% This function defines internal parameters ofor the robust control
% synthesis.

par.dt = par.T / (par.K - 1);  % Time step size

% Normalize initial and target states
par.X0=par.X0/norm(par.X0,2);          
par.XT=par.XT/norm(par.XT,2); 

% Legendre coefficients
par.alpha_under=(par.alpha_max-par.alpha_min)/2;
par.alpha_over=(par.alpha_max+par.alpha_min)/2;
par.beta_under=(par.beta_max-par.beta_min)/2;
par.beta_over=(par.beta_max+par.beta_min)/2;

% Legendre expansion coefficients
ka=0:1:par.Na-2;
ca=(ka+1)./((2*ka+3).*(2*ka+1)).^(1/2);
kb=0:1:par.Nb-2;
cb=(kb+1)./((2*kb+3).*(2*kb+1)).^(1/2);

% Expansion state and matrices
if (par.Na == 1 && par.Nb==1) || (par.alpha_under==0 && par.beta_under==0)
    par.A=par.A_base;
    par.B=par.B_base;
end

if par.alpha_under ~= 0 && par.beta_under ~= 0 && par.Na~=1 && par.Nb~=1
    par.X0 = [par.X0;spalloc(par.Na*par.Nb*par.n_base-par.n_base,1,0)];
    par.XT = [par.XT;spalloc(par.Na*par.Nb*par.n_base-par.n_base,1,0)];
    par.C_a = diag(par.alpha_over*ones(par.Na,1))+diag(par.alpha_under*ca,1)+diag(par.alpha_under*ca,-1); 
    par.C_b = diag(par.beta_over*ones(par.Nb,1))+diag(par.beta_under*cb,1)+diag(par.beta_under*cb,-1);
    par.C_a = sparse(par.C_a); par.C_b=sparse(par.C_b);
    par.A = kron(kron(par.C_a,speye(par.Nb)),par.A_base);
    par.B = [];
    for r = 1:par.m
        par.B = [par.B, kron(kron(speye(par.Na),par.C_b),par.B_base(:,(r-1)*par.n_base+1:r*par.n_base))];
    end
end

if (par.Na~=1 && par.Nb==1) || (par.alpha_under ~= 0 && par.beta_under == 0)
    par.X0 = [par.X0;spalloc(par.Na*par.n_base-par.n_base,1,0)];
    par.XT = [par.XT;spalloc(par.Na*par.n_base-par.n_base,1,0)];
    par.C_a = diag(par.alpha_over*ones(par.Na,1))+diag(par.alpha_under*ca,1)+diag(par.alpha_under*ca,-1); 
    par.C_a = sparse(par.C_a);
    par.A = kron(par.C_a,par.A_base);
    par.B = [];
    for r = 1:par.m
        par.B = [par.B, kron(speye(par.Na),par.B_base(:,(r-1)*par.n_base+1:r*par.n_base))];
    end
end

if (par.Na==1 && par.Nb~=1) || (par.alpha_under == 0 && par.beta_under ~= 0)
    par.X0 = [par.X0;spalloc((par.Nb)*par.n_base-par.n_base,1,0)];
    par.XT = [par.XT;spalloc((par.Nb)*par.n_base-par.n_base,1,0)];
    par.C_b = diag(par.beta_over*ones(par.Nb,1))+diag(par.beta_under*cb,1)+diag(par.beta_under*cb,-1);
    par.C_b = sparse(par.C_b);
    par.A = kron(speye(par.Nb),par.A_base);
    par.B = [];
    for r = 1:par.m
        par.B = [par.B, kron(par.C_b,par.B_base(:,(r-1)*par.n_base+1:r*par.n_base))];
    end
end

par.n=length(par.X0);

% Define projection matrix
if Schrodinger == true
    par.P=speye(par.n); par.P(par.n_target(end):par.n_base/2:end-par.n_base/2+par.n_target(end),:)=[];
else
    par.P=speye(par.n); par.P(par.n_target(end):par.n_base:end-par.n_base+par.n_target(end),:)=[];
end

% Finite difference matrix and vectorized control bounds
par.D=sparse(diag(ones(par.K-2,1),1)-diag(ones(par.K-1,1))); par.D(end,:)=[];           
par.D=kron(par.D,speye(par.m)); 
par.du_max=max(max(par.du_max))*ones(par.m*(par.K-2),1);                                              
par.du_min=min(min(abs(par.du_min)))*ones(par.m*(par.K-2),1);                                             
par.u_max=max(max(par.u_max))*ones(par.m*(par.K-1),1);                                              
par.u_min=min(min(abs(par.u_min)))*ones(par.m*(par.K-1),1);

% Initialize Jacobian
par.J = zeros(par.n,par.m*(par.K-1));


end