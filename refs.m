function ref = refs(wind_seg, p)
%DDP.REFS  Optimal operating trajectories for one wind segment.
%
%   ref = ddp.refs(wind_seg, p)
%
% Verbatim port of the reference block from Main_Control_Script.
%
%   Wr_opt = lambda_opt * v / R
%   Pa_opt = 0.5 * rho * pi * R^2 * v^3 * Cp_opt
%   Ta_opt = Pa_opt / Wr_opt
%   K_opt  = 0.5 * rho * pi * R^2 * R^3 * Cp_opt / lambda_opt^3
%   Tg_opt = K_opt * Wr_opt^2
%
% lambda_opt and Cp_opt are expanded to full-length row vectors because the
% rest of the solver indexes them per time step.

    wind_seg = reshape(wind_seg, 1, []);
    L = numel(wind_seg);

    ref.lambda_opt = p.lambda_opt * ones(1, L);
    ref.Cp_opt     = p.Cp_opt     * ones(1, L);

    ref.Wr_opt = (ref.lambda_opt .* wind_seg) / p.R;
    ref.Pa_opt = 0.5 * p.Rho * (pi * p.R^2) * wind_seg.^3 .* ref.Cp_opt;
    ref.Ta_opt = ref.Pa_opt ./ ref.Wr_opt;
    ref.K_opt  = 0.5 * p.Rho * (pi * p.R^2) * (p.R^3 .* ref.Cp_opt ./ ref.lambda_opt.^3);
    ref.Tg_opt = ref.K_opt .* ref.Wr_opt.^2;

    ref.Tg_max = max(ref.Tg_opt);
    ref.Tg_min = min(ref.Tg_opt);
    ref.L      = L;
end
