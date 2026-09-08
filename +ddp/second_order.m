function [f_xx, f_ux, f_uu] = second_order(x, wind, u)
%DDP.SECOND_ORDER  Second-order dynamics tensors for the DDP backward pass.
%   Verbatim port of second_order_dynamics from Main_Control_Script.
%
%   f_xx  4x4x4x4    f_ux  4x1x4    f_uu  4x1x1
    % Extract parameters
     R = 4.4;
     Kt = 52;
     Jt = 16;
     Rs = 0.01;
     Rho = 1.25;
     Ro = Rho;
     Ls = 0.01;
    % -------------------------- % % -------------------------- %
    % Initialize tensors
    f_xx = zeros(4,4,4,4);  % ∂²f_i/∂x_j∂x_k
    f_ux = zeros(4,1,4);    % ∂²f_i/∂u∂x_j
    f_uu = zeros(4,1,1);    % ∂²f_i/∂u²
    % -------------------------- % % -------------------------- %
    % State variables
    Wr = x(1);      % Rotor speed
    lambda = x(2);  % Tip-speed ratio
    Tg = x(3);      % Generator torque
    Cp = x(4);      % Power coefficient
    Beta = 0;       % Pitch angle (fixed for now)
    % -------------------------- % % -------------------------- %
    % Compute Cp derivatives (using your exact implementation)
    D = 0.08 * Beta + lambda;
    D = max(D, 0.01);
    expPart = exp((0.4375 / (Beta^3 + 1)) - (12.5 / D));
    % -------------------------- % % -------------------------- %
    % First derivative
    term1 = -25.52 / D^2;
    term2 = (25.52 / D - 0.088 * Beta - 0.1) * (12.5 / D^2);
    dCp_dlambda = expPart * (term1 + term2);
    % -------------------------- % % -------------------------- %
    % Second derivative
    piece01 = 25 - (156.25 / D);
    piece02 = 0.088*Beta + 1.1 + (0.8932 / (Beta^3 + 1)) - (25.52 / D);
    piece03 = 51.04 - (638 / D);
    piece04 = exp((0.4375 / (Beta^3 + 1)) - (12.5 / D));
    piece05 = D^3;
    d2Cp_dlambda2 = ((piece01 * piece02) + piece03) * piece04 / piece05;
    d2Cp_dlambda2 = min(max(d2Cp_dlambda2, -50), 50);
    % -------------------------- % % -------------------------- %
    % 1. f_xx terms (second derivatives of dynamics)
    
    % Rotor acceleration (Wr_dot) terms
    f_xx(1,1,1,1) = (3*Kt)/(Jt*Wr^2);  % From friction term
    % -------------------------- % % -------------------------- %
    % Tip-speed ratio (lambda_dot) terms
    f_xx(2,1,1,2) = R/(wind*Wr^2);
    f_xx(2,2,1,1) = R/(wind*Wr^2);
    % -------------------------- % % -------------------------- %
    % Cp_dot terms (most significant nonlinearities)
    coeff = (0.5*Rho*pi*R^2*wind^3)/(Jt*Wr^2);
    f_xx(4,2,2,2) = coeff * d2Cp_dlambda2 * (R/wind)^2;
    % -------------------------- % % -------------------------- %
    % Cross terms between Wr and lambda
    cross_term_Wr_lambda = coeff/Wr * (3*dCp_dlambda + lambda*d2Cp_dlambda2);
    f_xx(1,1,1,2) = cross_term_Wr_lambda;
    f_xx(1,1,2,1) = cross_term_Wr_lambda;
    % -------------------------- % % -------------------------- %
    % Generator torque coupling terms
    f_xx(1,3,1,3) = -1/(Jt*Ls);  % From Tg_dot dynamics
    % -------------------------- % % -------------------------- %
    % 2. f_ux terms (control-state mixed derivatives)
    
    % Rotor acceleration control coupling
    f_ux(1,1,1) = -1/(Jt*Ls);
    % -------------------------- % % -------------------------- %
    % Tip-speed ratio control coupling
    f_ux(2,1,1) = -R/(wind*Jt*Ls*Wr^2);
    % -------------------------- % % -------------------------- %
    % Cp dynamics control coupling
    f_ux(4,1,2) = -R/(Jt*Ls*wind) * dCp_dlambda;
    % -------------------------- % % -------------------------- %
    % Generator torque control coupling
    f_ux(3,1,3) = -Rs/(Ls^2);
    % -------------------------- % % -------------------------- %
    % 3. f_uu remains zero (control appears linearly)
    % -------------------------- % % -------------------------- %
    % Symmetry enforcement
    for i = 1:4
        for j = 1:4
            temp = squeeze(f_xx(i,j,:,:));
            f_xx(i,j,:,:) = 0.5*(temp + temp');
        end
    end
    % -------------------------- % % -------------------------- %
end
