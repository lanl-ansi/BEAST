function [x_bar,u_bar,par] = state_trasfer_control(par,Schrodinger,report_progress)

% This function attempts to compute the open-loop control action that 
% achieves the specified state transfer.

options = optimoptions('quadprog','Algorithm','interior-point-convex','LinearSolver','sparse','Display','none',...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-12);

% Initialize state and Jacobian
u_bar = par.u0(:,1:par.K-1);
x_bar = simulate(par,u_bar); 
par = Jacobian(x_bar,u_bar,par);

% Compute initial error
if Schrodinger
error = norm(par.P*(par.XT-x_bar(:,end)),2);
else
error = norm(par.XT-x_bar(:,end),2);
end
error_old=10*error;
error_opt = error; x_bar_opt = x_bar;  u_bar_opt = u_bar; 

% Sequantial quadratic program
 for iter = 1:par.iter_max
     u_bar = u_bar(:);
    
     % Compute control update
     if Schrodinger
     du = quadprog(2*sparse((par.P*par.J)'*(par.P*par.J))+2*error^2*par.lambda*par.dt*speye(par.m*(par.K-1)),2*(par.P*par.J)'*(par.P*(x_bar(1:par.n,end)-par.XT)),[par.D;-par.D;speye(par.m*(par.K-1));-speye(par.m*(par.K-1))],[par.du_max-par.D*u_bar;par.du_min+par.D*u_bar;par.u_max-u_bar;par.u_min+u_bar],[],[],[],[],0*u_bar,options);   
     else
     du = quadprog(2*sparse(par.J'*par.J)+2*error^2*par.lambda*par.dt*speye(par.m*(par.K-1)),2*par.J'*(x_bar(1:par.n,end)-par.XT),[par.D;-par.D;speye(par.m*(par.K-1));-speye(par.m*(par.K-1))],[par.du_max-par.D*u_bar;par.du_min+par.D*u_bar;par.u_max-u_bar;par.u_min+u_bar],[],[],[],[],0*u_bar,options);     
     end

     u_bar = u_bar + du;
     u_bar = reshape(u_bar,par.m,par.K-1);
     energy = .3*par.dt*u_bar(:)'*u_bar(:)+par.dt*u_bar(:)'*par.D'*par.D*u_bar(:);
    
     % Simulate state update
     x_old = x_bar;
     x_bar = simulate(par,u_bar);    
     if Schrodinger
     error = norm(par.P*(par.XT-x_bar(:,end)),2);
     else
     error = norm(par.XT-x_bar(:,end),2);
     end

     if report_progress == true
        disp(num2str(iter))
        fprintf('error=%d, energy=%d, max|du|=%d,  ', error, energy, max(abs(du(:))))
     end
    
     % Accept or reject update
     if error<=error_old 
        error_opt = error;
        x_bar_opt = x_bar;
        u_bar_opt = u_bar;

     % Calculate Jacobian 
        if mod(iter,2)==1
            par = Jacobian(x_bar,u_bar,par); % Second order accurate update
        else
            dx = x_bar(:,end)-x_old(:,end);
            par.J = par.J + (dx - par.J*du)*du'/(du'*du); % Broyden update
        end
     else        
         error = error_opt;
         x_bar = x_bar_opt;
         u_bar = u_bar_opt;
         par.lambda=2*par.lambda;
     end
     
     % Stopping criteria
     if (abs(error-error_old)<10^(-7) && error<error_old && error<.5) || error<par.epsilon || (max(abs(du(:)))<0.00001 && error<0.5)
        break
     end

   error_old = error;       
 end

end

