function Ta_dot = compute_Ta_dot(Wr, Wr_dot, v_wind, v_wind_dot, R, rho)
%DDP.COMPUTE_TA_DOT  Verbatim port of compute_Ta_dot from Main_Control_Script.
    % Current operating point
    lambda = R*Wr/max(v_wind, 0.1);  % Avoid division by zero
    Cp = 0.22*(116/lambda - 5)*exp(-12.5/lambda);  % Your Cp function
    
    % Partial derivatives
    dCp_dlambda = 0.22*( (-116/lambda^2)*exp(-12.5/lambda) + (116/lambda - 5)*(12.5/lambda^2)*exp(-12.5/lambda) );
    
    dlambda_dWr = R/v_wind;
    dlambda_dv_wind = -R*Wr/v_wind^2;
    
    % Chain rule
    Ta = 0.5*rho*pi*R^3 * (Cp/lambda) * v_wind^2;
    Ta_dot = Ta * ( (dCp_dlambda/Cp - 1/lambda)*(dlambda_dWr*Wr_dot + dlambda_dv_wind*v_wind_dot) + 2*v_wind_dot/v_wind );
end
