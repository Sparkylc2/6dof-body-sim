classdef Scene < handle
%SCENE 3d view of the rocket: the whole flight on the left, a close up on the right
%   s = bsim.rocket.Scene(P)                   P struct or packed
%   s = bsim.rocket.Scene(P, 'GimbalExaggeration', 3, 'Visible', 'on')
%   s.update(t, x, terms)                      x 14x1 state, terms 34x1 (or the struct)
%   s.update(t, x, terms, [], u)               u (5 + n_fins)x1 controls, so the fins get drawn
%                                              where they are (without u they're drawn stowed)
%   s.reset()                                  clears the trail
%
%   left:  the trail so far (orange while the engine's burning, blue after), the pad, where it is
%   right: the rocket around its cm (ENU, so up is up), with
%          orange plume   the exhaust leaving the nozzle, length goes with thrust
%          orange dashes  the thrust line, through the nozzle along the thrust direction, so the
%                         gap between it and the cm (black dot) is the moment arm
%          grey dashes    the body axis
%          blue arrow     the air going past (off to the side, dot end is downstream)
%          grey line      the last few seconds of the path
%          dark squares   the grid fins, folded flat against the body until deployed, then
%                         sticking out with the lattice across the flow, tilted by their deflection
%
%   GimbalExaggeration scales the drawn gimbal angle (only in the drawing) so small angles show
%   (the readout always gives the real ones), FinScale (default 3) draws the grid fins bigger
%   than they are. the dials top right show each fin's real deflection, seen from the side
%
%   bsim.rocket.replay and bsim.rocket.live are built on this

    properties
        Fig
        AxOverview
        AxClose
        GimbalExaggeration = 1
        TrailSeconds = 3 % how much path the close up shows
        FinScale = 3     % grid fins are drawn this much bigger than they are, so they show up
    end

    properties (Access = private)
        s_body   % body outline stations (forward from the dry cm) [m]
        r_body   % and radii [m]
        ring     % angles round the body
        L        % rocket length [m]
        T_ref    % thrust that draws a full length plume [N]
        x_G      % gimbal pivot station [m]
        L_n      % pivot to nozzle exit [m]
        R_body   % body radius [m]
        fin      % grid fin geometry (n, x, r, phi0, a = side of one fin)

        trail = zeros(0, 3) % ENU
        trail_t = zeros(0, 1)
        trail_burn = false(0, 1)

        h % graphics handles
    end

    methods
        function s = Scene(P, varargin)
            opts = inputParser;
            opts.addParameter('GimbalExaggeration', 1);
            opts.addParameter('Visible', 'on');
            opts.addParameter('Figure', []);
            opts.addParameter('FinScale', 3);
            opts.parse(varargin{:});
            s.GimbalExaggeration = opts.Results.GimbalExaggeration;
            s.FinScale = opts.Results.FinScale;

            if ~isstruct(P)
                P = bsim.rocket.unpack(P);
            end
            s.geometry(P);
            s.build(opts.Results.Figure, opts.Results.Visible);
        end

        function reset(s)
            s.trail = zeros(0, 3);
            s.trail_t = zeros(0, 1);
            s.trail_burn = false(0, 1);
        end

        function update(s, t, x, terms, history, u)
            %UPDATE draw the rocket at time t
            %   history (optional, [] to skip) is n x 5 [t pn pe pd T] and replaces the trail (for
            %   jumping around in a replay), otherwise this point is added to it
            %   u (optional) is the controls, for drawing the fins
            if isstruct(terms)
                terms = [terms.m; terms.x_cm; terms.Ix; terms.It; terms.T; terms.m_dot_e; ...
                    terms.V; terms.mach; terms.qbar; terms.alpha_t; terms.F_body; ...
                    terms.M_body; terms.F_aero; terms.M_aero; terms.g_i; terms.n_T; ...
                    terms.rho_e; terms.v_r];
            end
            x = x(:);
            terms = terms(:);
            if numel(terms) < 34
                error('bsim:error', ['terms has %d entries, the viewer needs the ones from ' ...
                    'bsim_rocket_xdot (refresh the xdot C Caller block so it picks up the ' ...
                    'new bsim_c.h, and restart matlab if the bridge was rebuilt while it was ' ...
                    'open)'], numel(terms));
            end
            if ~isvalid(s.Fig)
                return
            end

            pos = enu(x(1:3));
            burning = terms(5) > 0;
            if nargin < 6
                u = [];
            end
            if nargin >= 5 && ~isempty(history)
                s.trail = [history(:, 3), history(:, 2), -history(:, 4)];
                s.trail_t = history(:, 1);
                s.trail_burn = history(:, 5) > 0;
            elseif isempty(s.trail_t) || t > s.trail_t(end)
                s.trail(end + 1, :) = pos;
                s.trail_t(end + 1, 1) = t;
                s.trail_burn(end + 1, 1) = burning;
            end

            % --- overview ---
            b = s.trail_burn;
            on = b | [false; b(1:end - 1)]; % keeps the burnout point on both, so they join
            set(s.h.burn, 'XData', masked(s.trail(:, 1), on), 'YData', masked(s.trail(:, 2), on), ...
                'ZData', masked(s.trail(:, 3), on));
            set(s.h.coast, 'XData', masked(s.trail(:, 1), ~b), ...
                'YData', masked(s.trail(:, 2), ~b), 'ZData', masked(s.trail(:, 3), ~b));
            set(s.h.now, 'XData', pos(1), 'YData', pos(2), 'ZData', pos(3));

            % --- close up, everything relative to the cm in ENU ---
            R = quat_to_dcm(x(7:10)); % body -> NED
            B = ned_to_enu() * R;      % body -> ENU
            x_cm = terms(2);

            [S, Th] = ndgrid(s.s_body - x_cm, s.ring);
            Rr = repmat(s.r_body, 1, numel(s.ring));
            pts = B * [S(:)'; (Rr(:) .* cos(Th(:)))'; (Rr(:) .* sin(Th(:)))'];
            sz = size(S);
            set(s.h.body, 'XData', reshape(pts(1, :), sz), 'YData', reshape(pts(2, :), sz), ...
                'ZData', reshape(pts(3, :), sz));

            ax = B * [1; 0; 0];
            ends = [ax * (max(s.s_body) - x_cm + 0.15 * s.L), ax * (min(s.s_body) - x_cm - 0.15 * s.L)];
            set(s.h.axis, 'XData', ends(1, :), 'YData', ends(2, :), 'ZData', ends(3, :));

            n_T = terms(26:28);
            gp = atan2(n_T(3), n_T(1));
            gy = asin(max(-1, min(1, -n_T(2))));
            k = s.GimbalExaggeration;
            n_draw = [cos(k * gp) * cos(k * gy); -sin(k * gy); sin(k * gp) * cos(k * gy)];
            % the nozzle swings about the gimbal pivot, by the drawn angle (with no exaggeration
            % this lands exactly on rho_e, the model's nozzle exit, terms(29:31))
            pivot = [s.x_G - x_cm; 0; 0];
            exit_draw = pivot - s.L_n * n_draw;
            e = B * exit_draw;
            plume = B * (exit_draw - n_draw * (0.5 * s.L * terms(5) / s.T_ref));
            set(s.h.plume, 'XData', [e(1), plume(1)], 'YData', [e(2), plume(2)], ...
                'ZData', [e(3), plume(3)], 'Visible', onoff(burning));
            line_end = B * (exit_draw + n_draw * 1.2 * s.L);
            set(s.h.tline, 'XData', [e(1), line_end(1)], 'YData', [e(2), line_end(2)], ...
                'ZData', [e(3), line_end(3)]);
            set(s.h.nozzle, 'XData', e(1), 'YData', e(2), 'ZData', e(3));

            v_r = terms(32:34);
            if norm(v_r) > 1
                % drawn off to the side of the rocket (it'd be hidden inside the body otherwise),
                % pointing the way the air goes past, so from upstream towards the rocket
                side = cross(ax, [0; 0; 1]);
                if norm(side) < 0.1
                    side = cross(ax, [1; 0; 0]);
                end
                side = side / norm(side);
                d = B * (v_r / norm(v_r)) * 0.3 * s.L;
                tip = side * 0.25 * s.L;
                tail = tip + d;
                set(s.h.wind, 'XData', [tail(1), tip(1)], 'YData', [tail(2), tip(2)], ...
                    'ZData', [tail(3), tip(3)], 'Visible', 'on');
                set(s.h.windtip, 'XData', tip(1), 'YData', tip(2), 'ZData', tip(3), ...
                    'Visible', 'on');
                set(s.h.windlabel, 'Position', tail(:)', 'Visible', 'on');
            else
                set([s.h.wind, s.h.windtip, s.h.windlabel], 'Visible', 'off');
            end

            % --- grid fins ---
            deployed = false;
            delta = zeros(s.fin.n, 1);
            if numel(u) >= 5 + s.fin.n
                deployed = u(5) > 0.5;
                delta = u(6:5 + s.fin.n);
                delta = delta(:);
            end
            for k = 1:s.fin.n
                [corners, grid] = s.fin_shape(k, deployed, delta(k), x_cm);
                c = B * corners;
                g = B * grid;
                set(s.h.fin(k), 'XData', c(1, :)', 'YData', c(2, :)', 'ZData', c(3, :)');
                set(s.h.fingrid(k), 'XData', g(1, :), 'YData', g(2, :), 'ZData', g(3, :));

                % dial: the fin seen side on, turned by its deflection
                cx = 1.2 * (k - 1);
                d = 0.45 * [cos(delta(k)), sin(delta(k))];
                set(s.h.dial(k), 'XData', cx + [-d(1), d(1)], 'YData', [-d(2), d(2)]);
                if deployed
                    set(s.h.dial(k), 'Color', [0.85, 0.33, 0.1]);
                    s.h.dialtxt(k).String = sprintf('%+.1f', rad2deg(delta(k)));
                else
                    set(s.h.dial(k), 'Color', [0.6, 0.6, 0.6]);
                    s.h.dialtxt(k).String = 'stowed';
                end
            end

            recent = s.trail_t >= t - s.TrailSeconds & s.trail_t <= t;
            rel = s.trail(recent, :) - pos;
            set(s.h.recent, 'XData', rel(:, 1), 'YData', rel(:, 2), 'ZData', rel(:, 3));

            % --- readout ---
            if burning
                phase = 'burning';
            elseif -x(6) > 0
                phase = 'coasting up';
            else
                phase = 'coming down';
            end
            str = sprintf([ ...
                't  %7.2f s   (%s)\n' ...
                'h  %7.0f m     V %5.0f m/s   mach %4.2f\n' ...
                'climb %5.0f m/s   downrange %6.0f m\n' ...
                'T  %7.2f kN    m %6.1f kg\n' ...
                'gimbal  pitch %+5.2f deg   yaw %+5.2f deg\n' ...
                'alpha_t %6.2f deg   qbar %6.0f Pa'], ...
                t, phase, pos(3), terms(7), terms(8), -x(6), norm(pos(1:2)), terms(5) / 1e3, ...
                x(14), rad2deg(gp), rad2deg(gy), rad2deg(terms(10)), terms(9));
            if s.fin.n > 0
                if deployed
                    str = sprintf(['%s\nfins %s deg'], str, sprintf(' %+5.1f', rad2deg(delta)));
                else
                    str = sprintf('%s\nfins stowed', str);
                end
            end
            if s.GimbalExaggeration ~= 1
                str = sprintf('%s\n(gimbal drawn x%g)', str, s.GimbalExaggeration);
            end
            s.h.text.String = str;
        end
    end

    methods (Access = private)
        function [corners, grid] = fin_shape(s, k, deployed, delta, x_cm)
            % one grid fin in body axes relative to the cm: the square outline and a 3 x 3 lattice
            phi = s.fin.phi0 + 2 * pi * (k - 1) / s.fin.n;
            xb = [1; 0; 0];
            er = [0; cos(phi); sin(phi)]; % out from the body axis
            nk = [0; -sin(phi); cos(phi)]; % x cross er, across the fin
            a = s.fin.a * s.FinScale;
            if deployed
                % sticking out, lattice across the flow, turned about er by the deflection
                % (scaled up about its inner edge, so it still sits on the body)
                inner = s.fin.r - s.fin.a / 2;
                c = (s.fin.x - x_cm) * xb + (inner + a / 2) * er;
                d1 = er;
                d2 = cos(delta) * nk + sin(delta) * xb;
            else
                % folded back flat against the body, hinged at the fin station
                c = (s.fin.x - x_cm - a / 2) * xb + (s.R_body + 0.01) * er;
                d1 = xb;
                d2 = nk;
            end
            corners = c + d1 * ([-1, 1, 1, -1] * a / 2) + d2 * ([-1, -1, 1, 1] * a / 2);
            g = (-1:2 / 3:1) * a / 2;
            grid = zeros(3, 0);
            for v = g(2:end - 1)
                grid = [grid, c + d1 * v - d2 * a / 2, c + d1 * v + d2 * a / 2, nan(3, 1), ...
                    c + d2 * v - d1 * a / 2, c + d2 * v + d1 * a / 2, nan(3, 1)]; %#ok<AGROW>
            end
        end

        function geometry(s, P)
            sl = P.slices; % [x d dx]
            if isempty(sl)
                error('bsim:error', 'the rocket needs slices to draw its body');
            end
            [~, i] = sort(sl(:, 1));
            sl = sl(i, :);
            tail = sl(1, 1) - sl(1, 3) / 2;
            tip = sl(end, 1) + sl(end, 3) / 2;
            s.s_body = [tail; tail; sl(:, 1); tip];
            % closed at the base, and the nose comes to a point if the last slice is narrow
            tip_r = sl(end, 2) / 2;
            if tip_r < 0.2 * P.d
                tip_r = 0;
            end
            s.r_body = [0; sl(1, 2) / 2; sl(:, 2) / 2; tip_r];
            s.ring = linspace(0, 2 * pi, 20);
            s.L = tip - tail;
            s.x_G = P.engine.x_G;
            s.L_n = P.engine.L_n;
            s.R_body = P.d / 2;
            s.fin = struct('n', P.fins.n, 'x', P.fins.x, 'r', P.fins.r, 'phi0', P.fins.phi0, ...
                'a', sqrt(max(P.fins.S, 0)));

            if strcmp(P.engine.kind, 'thrust_curve') && ~isempty(P.engine.thrust_curve)
                s.T_ref = max(P.engine.thrust_curve(:, 2));
            else
                s.T_ref = P.engine.T_vac_max;
            end
            if ~(s.T_ref > 0)
                s.T_ref = 1;
            end
        end

        function build(s, fig, visible)
            if isempty(fig)
                fig = figure('Name', 'rocket', 'Color', 'w', 'Visible', visible, ...
                    'Position', [80, 80, 1300, 680]);
            end
            s.Fig = fig;
            orange = [0.93, 0.45, 0.1];
            blue = [0.15, 0.4, 0.8];

            % overview
            a = axes(fig, 'Position', [0.05, 0.14, 0.4, 0.78]);
            hold(a, 'on'); grid(a, 'on'); box(a, 'on'); view(a, 35, 20);
            % not to scale: the climb is km and the drift is often 100s of m, equal axes would
            % squash it into a line
            s.h.burn = plot3(a, nan, nan, nan, '-', 'Color', orange, 'LineWidth', 2);
            s.h.coast = plot3(a, nan, nan, nan, '-', 'Color', blue, 'LineWidth', 1.5);
            plot3(a, 0, 0, 0, 'k^', 'MarkerFaceColor', 'k');
            s.h.now = plot3(a, 0, 0, 0, 'o', 'Color', 'k', 'MarkerFaceColor', orange, ...
                'MarkerSize', 7);
            xlabel(a, 'east [m]'); ylabel(a, 'north [m]'); zlabel(a, 'up [m]');
            title(a, 'flight');
            legend(a, {'powered', 'coast', 'pad'}, 'Location', 'northeast', 'AutoUpdate', 'off');
            s.AxOverview = a;

            % close up
            c = axes(fig, 'Position', [0.52, 0.14, 0.45, 0.72]);
            hold(c, 'on'); grid(c, 'on'); box(c, 'on'); view(c, 35, 15);
            lim = 0.75 * s.L;
            axis(c, [-lim, lim, -lim, lim, -lim, lim]);
            c.XLimMode = 'manual'; c.YLimMode = 'manual'; c.ZLimMode = 'manual';
            daspect(c, [1, 1, 1]); pbaspect(c, [1, 1, 1]); % to scale, and stays put
            s.h.recent = plot3(c, nan, nan, nan, '-', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 1);
            s.h.body = surf(c, nan(2), nan(2), nan(2), 'FaceColor', [0.82, 0.83, 0.86], ...
                'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.5);
            s.h.axis = plot3(c, nan, nan, nan, '--', 'Color', [0.45, 0.45, 0.45]);
            s.h.tline = plot3(c, nan, nan, nan, '--', 'Color', orange, 'LineWidth', 1.2);
            s.h.plume = plot3(c, nan, nan, nan, '-', 'Color', orange, 'LineWidth', 6);
            s.h.nozzle = plot3(c, nan, nan, nan, 's', 'Color', 'k', 'MarkerFaceColor', orange);
            s.h.wind = plot3(c, nan, nan, nan, '-', 'Color', blue, 'LineWidth', 2);
            s.h.windtip = plot3(c, nan, nan, nan, 'o', 'Color', blue, 'MarkerFaceColor', blue);
            s.h.windlabel = text(c, 0, 0, 0, ' air', 'Color', blue, 'Visible', 'off');
            plot3(c, 0, 0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            s.h.fin = gobjects(s.fin.n, 1);
            s.h.fingrid = gobjects(s.fin.n, 1);
            for k = 1:s.fin.n
                s.h.fin(k) = patch(c, nan(4, 1), nan(4, 1), nan(4, 1), [0.25, 0.27, 0.32], ...
                    'FaceAlpha', 0.35, 'EdgeColor', [0.15, 0.15, 0.18], 'LineWidth', 1.2);
                s.h.fingrid(k) = plot3(c, nan, nan, nan, '-', 'Color', [0.15, 0.15, 0.18]);
            end
            light(c, 'Position', [1, -1, 2]);
            xlabel(c, 'east [m]'); ylabel(c, 'north [m]'); zlabel(c, 'up [m]');
            title(c, 'close up (around the cm)');
            s.h.text = text(c, 0.02, 0.98, '', 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'FontName', 'Menlo', 'FontSize', 10, 'BackgroundColor', [1, 1, 1, 0.8]);
            s.AxClose = c;

            % fin dials, top right above the close up
            d = axes(fig, 'Position', [0.8, 0.86, 0.18, 0.11], 'Visible', 'off');
            hold(d, 'on'); axis(d, 'equal');
            axis(d, [-0.7, 1.2 * max(s.fin.n - 1, 0) + 0.7, -0.75, 0.6]);
            s.h.dial = gobjects(s.fin.n, 1);
            s.h.dialtxt = gobjects(s.fin.n, 1);
            for k = 1:s.fin.n
                cx = 1.2 * (k - 1);
                plot(d, cx + [-0.5, 0.5], [0, 0], ':', 'Color', [0.7, 0.7, 0.7]);
                s.h.dial(k) = plot(d, cx + [-0.45, 0.45], [0, 0], '-', 'LineWidth', 4, ...
                    'Color', [0.6, 0.6, 0.6]);
                text(d, cx, 0.45, sprintf('fin %d', k), 'HorizontalAlignment', 'center', ...
                    'FontSize', 8, 'Color', [0.3, 0.3, 0.3]);
                s.h.dialtxt(k) = text(d, cx, -0.6, '', 'HorizontalAlignment', 'center', ...
                    'FontSize', 8, 'FontName', 'Menlo');
            end
        end
    end
end

function v = enu(ned)
v = [ned(2), ned(1), -ned(3)];
end

function M = ned_to_enu()
M = [0, 1, 0; 1, 0, 0; 0, 0, -1];
end

function R = quat_to_dcm(q)
% body -> NED, from the [e0 e1 e2 e3] attitude (same as frames.cpp)
q = q / norm(q);
e0 = q(1); e1 = q(2); e2 = q(3); e3 = q(4);
R = [e1^2 + e0^2 - e2^2 - e3^2, 2 * (e1 * e2 - e3 * e0), 2 * (e1 * e3 + e2 * e0);
    2 * (e1 * e2 + e3 * e0), e2^2 + e0^2 - e1^2 - e3^2, 2 * (e2 * e3 - e1 * e0);
    2 * (e1 * e3 - e2 * e0), 2 * (e2 * e3 + e1 * e0), e3^2 + e0^2 - e1^2 - e2^2];
end

function v = masked(v, keep)
v(~keep) = nan;
end

function s = onoff(b)
if b, s = 'on'; else, s = 'off'; end
end
