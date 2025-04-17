function [x_bar,u_bar,x,par] = solver(par,initialize,Schrodinger,min_energy,plot_sol,report_progress)

% This function organizes the operations of the primary functions.

% Initialize control
if par.Na==1 && par.Nb==1
    initialize = false;
end

if initialize
    Na = par.Na;
    Nb = par.Nb;
    par.Na = 1;
    par.Nb = 1;
    par = param(par,Schrodinger);
    [x_bar,u_bar,par] = state_trasfer_control(par,Schrodinger,report_progress);
    [~,u_bar,~] = min_energy_control(par,x_bar,u_bar,report_progress);
    par.Na = Na;
    par.Nb = Nb;
    par = param(par,Schrodinger);
    par.u0 = u_bar;
end

if initialize == false
    par = param(par,Schrodinger);
end

% State transfer control
[x_bar,u_bar,par] = state_trasfer_control(par,Schrodinger,report_progress);
% x = x_bar;

% Minimal energy control
if min_energy == true
    [x_bar,u_bar,par] = min_energy_control(par,x_bar,u_bar,report_progress);
end

% Figures
if plot_sol == true
    x = plots(par,x_bar,u_bar,Schrodinger);
end

end

