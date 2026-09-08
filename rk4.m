function x_next = rk4(x, u, Tg_des_k, wind_k, k, dt)
%DDP.RK4  One RK4 step of the 4-state turbine model. Verbatim port of
%   turbine_dynamics_rk4 from Main_Control_Script.
%
%   x = [Wr; lambda; Cp; Tg],  u scalar,  returns x_next [4x1]
%
% PRESERVED AS-IS, including:
%   * plant constants hardcoded here rather than taken from ddp.params
%     (the original did the same; changing it would risk divergence)
%   * Cp uses  0.22*(116/lambda - 0.4*beta - 5)*exp(-12.5/lambda)
%     in the RK4 stages, while dCp/dlambda uses the OTHER model,
%     exp(0.4375/(beta^3+1) - 12.5/D). The two disagree; see README.
%   * the soft coupling  Tg_next += 0.25*(Tg_des_k - Tg_next)
     R = 4.4;
     Kt = 52;
     Jt = 16;
     Rs = 0.01;
     Rho = 1.25;
     Ro = Rho;
     Ls = 0.01;
    % -------------------------- % % -------------------------- %
    % x:     [Wr; lambda; Cp; Tg]
    % u:     scalar input
    % wind_k: current wind value
    % dt:    time step
    % Returns: x_next [4×1] after 1 RK4 step
    % -------------------------- % % -------------------------- %
    wind_k = max(wind_k, 0.1);
    % -------------------------- % % -------------------------- %
    % --- Unpack state ---
    Wr_k     = x(1);
    lambda_k = x(2);
    lambda_k = max(real(lambda_k), 0.5);
    Cp_k     = x(3);
    Tg_k     = x(4);
    Beta     = 0;
    tau_d = 0.1 * sin(2 * pi * k * dt);
    U_k = u;
    % -------------------------- % % -------------------------- %
    % --- Disturbance (optional) ---
    % tau_d = 0;  % set to zero or pass as argument if needed
    % -------------------------- % % -------------------------- %
    %                   RK4 state Updating 
    % -------------------------- % % -------------------------- %
    % compute dCp/dlambda at current state
    Beta = 0;
    D = max(0.08*Beta + lambda_k, 1e-2);
    expPart   = exp((0.4375/(Beta^3+1)) - (12.5/D));
    term1     = -25.52 / D^2;
    term2     = (25.52/D - 0.088*Beta - 0.1)*(12.5/D^2);
    dCp_dlam1 = expPart*(term1 + term2);
    % -------------------------- % % -------------------------- %
    % --- RK4 Step 1 ---
    Cp1    = 0.22*((116/lambda_k) - 0.4*Beta - 5)*exp(-12.5/lambda_k);
    Ta1    = 0.5*Ro*pi*R^3*wind_k^2*(Cp1/lambda_k);
    Wr_dot1    = (1/Jt)*(Ta1 - Kt*Wr_k + (Rs/Ls)*Tg_k - (1/Ls)*U_k + tau_d);
    lambda_dot1 = (R/wind_k)*Wr_dot1;
    Cp_dot1     = dCp_dlam1 * lambda_dot1;
    Tg_dot1     = (1/Ls)*(-Rs*Tg_k + U_k);
    % -------------------------- % % -------------------------- %
    % --- RK4 Step 2 ---
    Wr2     = Wr_k     + 0.5*dt*Wr_dot1;
    lambda2 = (R*Wr2)/wind_k;
    Cp2     = 0.22*((116/lambda2) - 0.4*Beta - 5)*exp(-12.5/lambda2);
    Tg2     = Tg_k     + 0.5*dt*Tg_dot1;
    tau_d2  = 0.1 * sin(2*pi*(k+0.5)*dt);
    % recalc derivative
    D2        = max(0.08*Beta + lambda2, 1e-2);
    expP2     = exp((0.4375/(Beta^3+1)) - (12.5/D2));
    term12    = -25.52/D2^2;
    term22    = (25.52/D2 -0.088*Beta -0.1)*(12.5/D2^2);
    dCp_dlam2 = expP2*(term12+term22);
    % -------------------------- % % -------------------------- %
    Ta2        = 0.5*Ro*pi*R^3*wind_k^2*(Cp2/lambda2);
    Wr_dot2    = (1/Jt)*(Ta2 - Kt*Wr2 + (Rs/Ls)*Tg2 - (1/Ls)*U_k + tau_d2);
    lambda_dot2 = (R/wind_k)*Wr_dot2;
    Cp_dot2     = dCp_dlam2 * lambda_dot2;
    Tg_dot2     = (1/Ls)*(-Rs*Tg2 + U_k);
    % -------------------------- % % -------------------------- %
    % --- RK4 Step 3 ---
    Wr3     = Wr_k     + 0.5*dt*Wr_dot2;
    lambda3 = (R*Wr3)/wind_k;
    Cp3     = 0.22*((116/lambda3) - 0.4*Beta - 5)*exp(-12.5/lambda3);
    Tg3     = Tg_k     + 0.5*dt*Tg_dot2;
    tau_d3  = tau_d2;  % same half‐step disturbance
    % -------------------------- % % -------------------------- %
    % recalc derivative
    D3        = max(0.08*Beta + lambda3, 1e-2);
    expP3     = exp((0.4375/(Beta^3+1)) - (12.5/D3));
    term13    = -25.52/D3^2;
    term23    = (25.52/D3 -0.088*Beta -0.1)*(12.5/D3^2);
    dCp_dlam3 = expP3*(term13+term23);
    % -------------------------- % % -------------------------- %
    Ta3        = 0.5*Ro*pi*R^3*wind_k^2*(Cp3/lambda3);
    Wr_dot3    = (1/Jt)*(Ta3 - Kt*Wr3 + (Rs/Ls)*Tg3 - (1/Ls)*U_k + tau_d3);
    lambda_dot3 = (R/wind_k)*Wr_dot3;
    Cp_dot3     = dCp_dlam3 * lambda_dot3;
    Tg_dot3     = (1/Ls)*(-Rs*Tg3 + U_k);
    % -------------------------- % % -------------------------- %
    % --- RK4 Step 4 ---
    Wr4     = Wr_k     + dt*Wr_dot3;
    lambda4 = (R*Wr4)/wind_k;
    Cp4     = 0.22*((116/lambda4) - 0.4*Beta - 5)*exp(-12.5/lambda4);
    Tg4     = Tg_k     + dt*Tg_dot3;
    tau_d4  = 0.1 * sin(2*pi*(k+1)*dt);
    % -------------------------- % % -------------------------- %
    % recalc derivative
    D4        = max(0.08*Beta + lambda4, 1e-2);
    expP4     = exp((0.4375/(Beta^3+1)) - (12.5/D4));
    term14    = -25.52/D4^2;
    term24    = (25.52/D4 -0.088*Beta -0.1)*(12.5/D4^2);
    dCp_dlam4 = expP4*(term14+term24);
    % -------------------------- % % -------------------------- %
    Ta4        = 0.5*Ro*pi*R^3*wind_k^2*(Cp4/lambda4);
    Wr_dot4    = (1/Jt)*(Ta4 - Kt*Wr4 + (Rs/Ls)*Tg4 - (1/Ls)*U_k + tau_d4);
    lambda_dot4 = (R/wind_k)*Wr_dot4;
    Cp_dot4     = dCp_dlam4 * lambda_dot4;
    Tg_dot4     = (1/Ls)*(-Rs*Tg4 + U_k);
    % -------------------------- % % -------------------------- %
    % State update
    Tg_next = Tg_k + (dt/6)*(Tg_dot1 + 2*Tg_dot2 + 2*Tg_dot3 + Tg_dot4);
    %     if mod(k,10) == 0  % Only update Wr every 10 steps
    %         Wr_next = Wr_k + (dt/6)*(Wr_dot1 + 2*Wr_dot2 + 2*Wr_dot3 + Wr_dot4);
    %     else
    %         % Wr_next = Wr_k;
    %         Ta = (0.5 * Ro * pi * R^3) * (wind_k^2) * (Cp_k / lambda_k);
    %         Wr_next = Wr_k + dt*(Ta - Kt*Wr_k - Tg_k + tau_d)/Jt; 
    %     end
    Wr_next = Wr_k + (dt/6)*(Wr_dot1 + 2*Wr_dot2 + 2*Wr_dot3 + Wr_dot4);
    
    %Tg_next = Tg_k + (1/Ls)*(-Rs*Tg_k + U_k);
    % --- MODIFIED --- Soft constraint coupling ( necessary for numerical
    % stability)
    Tg_next = Tg_next + (0.25)*(Tg_des_k - Tg_next);
    % Recompute derived states
    lambda_k1 = (R*Wr_next)/wind_k;
    Cp_k1 = 0.22*((116/lambda_k1) - 5)*exp(-12.5/lambda_k1);
    % save into X
    % x_next = [Wr_next; lambda_k1; Cp_k1; Tg_next];
    % -------------------------- % % -------------------------- %
    % Cp_k1 = max( min( Cp_k1 , 0.5) , 0 );
    % lambda_k1 = max( min( lambda_k1 , 15) , 2 );
    % Tg_next = max( min( Tg_next , 2e3) , 0 );
    x_next = [Wr_next; lambda_k1; Cp_k1; Tg_next];
    % -------------------------- % % -------------------------- %
end

%%
% ------------------------ % % ------------------------ %
% Function #1 wind turbine dynamics + RK4 state updating
% ------------------------ % % ------------------------ %
