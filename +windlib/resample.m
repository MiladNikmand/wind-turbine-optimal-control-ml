function [x_out, y_out] = resample(y_in, step_size, T_end)
%WINDLIB.RESAMPLE  Re-grid a raw wind trace onto a uniform step_size grid.
%
%   [x, y] = windlib.resample(y_in, step_size, T_end)
%
% Faithful port of Intp_data(..., start_condition = 1) from
% wind_profile_loader.m. Behaviour is unchanged:
%
%   * the original x column of the digitised data is DISCARDED -- time is
%     rebuilt as linspace(0, T_end, numel(y_in)). The digitiser x values in
%     the source CSVs are plot-trace artefacts and were never used.
%   * pchip interpolation (shape preserving, no overshoot at the knots)
%   * desired_points = T_end/step_size + 1, so the returned grid spacing is
%     exactly step_size whenever T_end/step_size is an integer.
%
% Returns row vectors, matching interp1's behaviour when the query grid is
% a row -- the rest of project 1 assumes `wind` is a row vector.

    y_in = y_in(:);

    x_original = linspace(0, T_end, numel(y_in));
    y_original = y_in;

    desired_points = (T_end ./ step_size) + 1;

    if abs(desired_points - round(desired_points)) > 1e-9
        warning('windlib:resample:nonIntegerGrid', ...
            ['T_end/step_size = %.6f is not an integer, so the resampled grid ' ...
             'spacing will be %.8g instead of the requested %.8g.'], ...
            T_end/step_size, T_end/(round(desired_points)-1), step_size);
    end
    desired_points = round(desired_points);

    x_out = linspace(min(x_original), max(x_original), desired_points);
    y_out = interp1(x_original, y_original, x_out, 'pchip');

    x_out = reshape(x_out, 1, []);
    y_out = reshape(y_out, 1, []);
end
