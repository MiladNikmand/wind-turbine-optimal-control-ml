function [f_xx, f_ux, f_uu] = second_order_v2(x, wind, u, p)
%DDP.SECOND_ORDER_V2  Second-order dynamics tensors, derived from scratch.
%
%   [f_xx, f_ux, f_uu] = ddp.second_order_v2(x, wind, u, p)
%
% Drop-in replacement for ddp.second_order. Same output shapes, so
% ddp.backward consumes it unchanged. Two differences:
%
%   1. STATE ORDERING.  Uses x = [Wr; lambda; Cp; Tg], matching ddp.rk4,
%      ddp.forward, ddp.backward and ddp.params. The original reads
%      x(3) as Tg and x(4) as Cp, and writes its tensor rows to match
%      that reversed convention -- see test_second_order.
%
%   2. NO HARDCODED CONSTANTS.  Takes p from ddp.params, so editing the
%      rotor radius in one place actually changes the backward pass.
%      p is optional; it defaults to ddp.params().
%
% DERIVATION.  With  K = 0.5*Rho*pi*R^3  and  g(lambda) = dCp/dlambda:
%
%   f1 = Wr_dot     = (1/Jt)*( K*wind^2*Cp/lambda - Kt*Wr + (Rs/Ls)*Tg - U/Ls )
%   f2 = lambda_dot = (R/wind) * f1
%   f3 = Cp_dot     = g(lambda) * f2
%   f4 = Tg_dot     = (1/Ls)*( -Rs*Tg + U )
%
% f1 is nonlinear only through Cp/lambda. f2 inherits that, scaled. f3 adds
% the product rule against g. f4 is affine, so its whole row is zero -- the
% property test_second_order checks first.
%
% The third derivative d3Cp/dlambda3, needed for one term of f3, has no
% clean closed form here and is taken by central difference of the analytic
% second derivative. Its contribution is small at sensible lambda.
%
% CAUTION: this reflects the aerodynamic torque used in ddp.rk4,
%   Ta = 0.5*Rho*pi*R^3*wind^2*Cp/lambda.
% ddp.backward's analytic A matrix implies Ta with lambda^3 in the
% denominator instead. Those two disagree independently of the ordering
% question. Decide which torque you intend before trusting either.

    if nargin < 4 || isempty(p)
        p = ddp.params();
    end

    R   = p.R;
    Kt  = p.Kt;
    Jt  = p.Jt;
    Rs  = p.Rs;
    Rho = p.Rho;
    Ls  = p.Ls;

    f_xx = zeros(4,4,4,4);
    f_ux = zeros(4,1,4);
    f_uu = zeros(4,1,1);

    Wr     = x(1);
    lambda = max(real(x(2)), 0.5);
    Cp     = x(3);
    %Tg    = x(4);   % affine throughout; never needed for curvature

    K  = 0.5 * Rho * pi * R^3;
    c  = K * wind^2 / Jt;      % coefficient on the Cp/lambda term of f1
    sc = R / wind;             % f2 = sc * f1

    % ---- Cp derivatives w.r.t. lambda ----
    g1 = dCp_dlambda(lambda);          % dCp/dlambda
    g2 = d2Cp_dlambda2(lambda);        % d2Cp/dlambda2
    hstep = 1e-4 * max(1, abs(lambda));
    g3 = (d2Cp_dlambda2(lambda + hstep) - d2Cp_dlambda2(lambda - hstep)) / (2*hstep);

    % =====================================================================
    %  f1 = Wr_dot
    % =====================================================================
    H1 = zeros(4,4);
    H1(2,2) =  2 * c * Cp / lambda^3;     % d2/dlambda2
    H1(2,3) = -c / lambda^2;              % d2/dlambda dCp
    H1(3,2) = H1(2,3);

    % =====================================================================
    %  f2 = (R/wind) * f1
    % =====================================================================
    H2 = sc * H1;

    % first derivatives of f2, needed by the product rule below
    df2 = zeros(4,1);
    df2(1) = sc * (-Kt / Jt);
    df2(2) = sc * (-c * Cp / lambda^2);
    df2(3) = sc * ( c / lambda);
    df2(4) = sc * ( Rs / (Jt*Ls));
    f2_val = sc * (1/Jt) * (K*wind^2*Cp/lambda - Kt*Wr + (Rs/Ls)*x(4) - u/Ls);

    % =====================================================================
    %  f3 = g(lambda) * f2      (product rule)
    % =====================================================================
    H3 = g1 * H2;
    H3(2,2) = H3(2,2) + g3 * f2_val + 2 * g2 * df2(2);
    for k = [1 3 4]
        H3(2,k) = H3(2,k) + g2 * df2(k);
        H3(k,2) = H3(2,k);
    end

    % =====================================================================
    %  f4 = Tg_dot  -- affine, so H4 stays zero
    % =====================================================================
    H4 = zeros(4,4);

    % ---- pack. ddp.backward sums f_xx(i,j,:,:) over j, so each Hessian
    %      goes in a single j slot and the rest stay zero. ----
    f_xx(1,1,:,:) = H1;
    f_xx(2,1,:,:) = H2;
    f_xx(3,1,:,:) = H3;
    f_xx(4,1,:,:) = H4;

    % =====================================================================
    %  mixed control-state derivatives
    % =====================================================================
    % f1 and f2 are linear in u, so d2/du dx = 0 for both.
    % f3 = g(lambda)*f2 picks up one term: d2f3/du dlambda = g2 * df2/du.
    df2_du = sc * (-1/(Jt*Ls));
    f_ux(3,1,2) = g2 * df2_du;
    % f4 linear in u and x -> zero.

    % f_uu is zero: the control enters every equation linearly.

    % ---- symmetry, matching the original's final step ----
    for i = 1:4
        for j = 1:4
            tmp = squeeze(f_xx(i,j,:,:));
            f_xx(i,j,:,:) = 0.5 * (tmp + tmp.');
        end
    end
end

% =========================================================================
function g = dCp_dlambda(lambda)
    Beta = 0;
    D    = max(0.08*Beta + lambda, 1e-2);
    expP = exp((0.4375/(Beta^3+1)) - (12.5/D));
    g    = expP * (-25.52/D^2 + (25.52/D - 0.088*Beta - 0.1)*(12.5/D^2));
end

function g2 = d2Cp_dlambda2(lambda)
    Beta = 0;
    D    = max(0.08*Beta + lambda, 1e-2);
    p01  = 25 - (156.25 / D);
    p02  = 0.088*Beta + 1.1 + (0.8932/(Beta^3+1)) - (25.52/D);
    p03  = 51.04 - (638 / D);
    p04  = exp((0.4375/(Beta^3+1)) - (12.5/D));
    g2   = ((p01 * p02) + p03) * p04 / D^3;
    g2   = min(max(g2, -50), 50);      % same clamp as the original
end