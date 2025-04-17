function [x_bar,u_bar,par] = min_energy_control(par,x_bar,u_bar,report_progress)

% This function iteratively updates the control action computed with the
% state transfer control to achieve a control action of minimal energy that
% performs the same state transfer. This function is bypassed if
% min_energy = false.

 options = optimoptions('quadprog','Algorithm','interior-point-convex','LinearSolver','sparse', ...
           'Display','none','OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12);

 % Energy of state transfer control
 energy = .3*par.dt*u_bar(:)'*u_bar(:)+par.dt*u_bar(:)'*par.D'*par.D*u_bar(:);
 energy_old = energy;
 XT = x_bar(:,end);

 % Sequantial quadratic program
 for  iter = 1:par.iter_max
        u_bar = u_bar(:);
        x_old = x_bar;

        % Control update
        du = quadprog(par.dt*(.09+par.mu)*speye(par.m*(par.K-1))+par.dt*par.D'*par.D, par.dt*(.3*speye(par.m*(par.K-1))+par.D'*par.D)*u_bar, [-speye(par.m*(par.K-1));speye(par.m*(par.K-1))], [par.u_min+u_bar;par.u_max-u_bar], par.P*par.J,par.P*(XT-x_bar(:,end)),[],[],0*u_bar,options);     
        energy = .3*par.dt*(u_bar(:)+du)'*(u_bar(:)+du)+par.dt*(u_bar(:)+du)'*par.D'*par.D*(u_bar(:)+du);
        
        % Stopping criteria
        if energy>energy_old || max(abs(du(:)))>0.5 || max(abs(du(:)))<0.001
            x_bar = x_old;
            u_bar = reshape(u_bar,par.m,par.K-1);
            break
        end
        
        % Control and corresponding state update
        u_bar = reshape(u_bar + du,par.m,par.K-1);
        x_old = x_bar;
        x_bar = simulate(par,u_bar);
        error = norm(par.P*(par.XT-x_bar(:,end)),2);  
        energy = .3*par.dt*u_bar(:)'*u_bar(:)+par.dt*u_bar(:)'*par.D'*par.D*u_bar(:);

        if report_progress == true
            disp(num2str(iter)) 
            fprintf('error=%d, energy=%d, max|du|=%d,  ', error, energy, max(abs(du(:))))
        end
        
        % Jacobian update
        if mod(iter,2)==1
            par = Jacobian(x_bar,u_bar,par); % Second order accurate update
        else
            dx = x_bar(:,end)-x_old(:,end);
            par.J = par.J + (dx - par.J*du)*du'/(du'*du); % Broyden update
        end
        energy_old = energy;
 end
 
end

