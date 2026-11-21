function varargout = replay(varargin)
%REPLAY play back a logged flight in the 3d view, with a time slider
%   bsim.rocket.replay(out, P)                 out from sim(), logged as out.xlog, out.terms
%                                              and (if it's there) out.ulog for the fins
%   bsim.rocket.replay(out, P, 'X', 'xlog', 'Terms', 'terms', 'U', 'ulog')   other names
%   bsim.rocket.replay(xlog, terms, P)         the logged timeseries themselves
%   bsim.rocket.replay(xlog, terms, P, 'U', ulog)
%   bsim.rocket.replay(t, X, TV, P)            t n x 1, X n x 14, TV n x 34
%   bsim.rocket.replay(t, X, TV, P, 'U', Um)   Um n x (5 + fins), the controls
%   ... 'Speed', 2                             play at 2x (default 1)
%   ... 'GimbalExaggeration', 3                draw the gimbal angles 3x bigger
%   ... 'Autoplay', false                      open paused at the start
%   scene = bsim.rocket.replay(...)            the bsim.rocket.Scene, to drive it yourself
%
%   play/pause with the button, drag the slider to look at any moment, the readout in the close
%   up shows the real numbers at that time (see bsim.rocket.Scene for what's drawn)
%   it returns straight away and plays in the background (matlab stays usable), so it's fine to
%   call from a callback like the model's StopFcn
Um = [];
if isa(varargin{1}, 'Simulink.SimulationOutput')
    out = varargin{1};
    P = varargin{2};
    opts = parse(varargin(3:end), 'ulog');
    xs = out.get(opts.X);
    ts = out.get(opts.Terms);
    us = [];
    if ischar(opts.U) && any(strcmp(out.who, opts.U))
        us = out.get(opts.U);
    end
    [t, X, TV, Um, Ps] = from_timeseries(xs, ts, us, P);
elseif isa(varargin{1}, 'timeseries')
    [xs, ts, P] = varargin{1:3};
    opts = parse(varargin(4:end), []);
    [t, X, TV, Um, Ps] = from_timeseries(xs, ts, opts.U, P);
else
    [t, X, TV, P] = varargin{1:4};
    opts = parse(varargin(5:end), []);
    t = t(:);
    Um = opts.U;
    Ps = P;
end
if isempty(Um)
    Um = zeros(numel(t), 0);
end

scene = bsim.rocket.Scene(Ps, 'GimbalExaggeration', opts.GimbalExaggeration);
fig = scene.Fig;
hist = [t, X(:, 1:3), TV(:, 5)];

% everything the controls need, kept on the figure
st.k = 1;
st.playing = false;
fig.UserData = st;

slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.12, 0.025, 0.7, 0.035], 'Min', t(1), 'Max', t(end), 'Value', t(1), ...
    'SliderStep', [0.002, 0.05] * max(1, (t(end) - t(1)) / 100));
button = uicontrol(fig, 'Style', 'pushbutton', 'Units', 'normalized', ...
    'Position', [0.03, 0.02, 0.07, 0.045], 'String', 'play');
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [0.84, 0.022, 0.05, 0.035], ...
    'String', 'speed', 'BackgroundColor', 'w');
speed = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', ...
    'Position', [0.89, 0.022, 0.07, 0.04], 'String', {'0.25x', '0.5x', '1x', '2x', '5x', '10x'}, ...
    'Value', 3);
speeds = [0.25, 0.5, 1, 2, 5, 10];
[~, i] = min(abs(speeds - opts.Speed));
speed.Value = i;

    function show(k)
        k = max(1, min(numel(t), k));
        st = fig.UserData;
        st.k = k;
        fig.UserData = st;
        scene.update(t(k), X(k, :), TV(k, :), hist(1:k, :), Um(k, :));
        slider.Value = t(k);
    end

% playback runs off a repeating timer, one frame per tick, so matlab stays free while it plays
t0 = t(1);
clock0 = tic;
player = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.033, 'BusyMode', 'drop', ...
    'TimerFcn', @(~, ~) tick());
fig.DeleteFcn = @(~, ~) cleanup();

    function play()
        st = fig.UserData;
        if st.k >= numel(t)
            show(1);
            st = fig.UserData;
        end
        st.playing = true;
        fig.UserData = st;
        button.String = 'pause';
        t0 = t(st.k);
        clock0 = tic;
        if strcmp(player.Running, 'off')
            start(player);
        end
    end

    function stop_playing()
        if isvalid(player)
            stop(player);
        end
        if isvalid(fig)
            st = fig.UserData;
            st.playing = false;
            fig.UserData = st;
            button.String = 'play';
        end
    end

    function tick()
        try
            if ~isvalid(fig) || ~fig.UserData.playing
                stop_playing();
                return
            end
            tn = t0 + toc(clock0) * speeds(speed.Value);
            k = find(t <= tn, 1, 'last');
            show(k);
            drawnow limitrate
            if k >= numel(t)
                stop_playing();
            end
        catch err
            stop_playing();
            if isvalid(fig)
                warning('bsim:replay', 'replay stopped: %s', err.message);
            end
        end
    end

    function cleanup()
        if isvalid(player)
            stop(player);
            delete(player);
        end
    end

slider.Callback = @(~, ~) slide();
button.Callback = @(~, ~) toggle();

    function slide()
        stop_playing();
        show(find(t <= slider.Value, 1, 'last'));
    end

    function toggle()
        if fig.UserData.playing
            stop_playing();
        else
            play();
        end
    end

% fit the overview round the whole flight so it doesn't jump about while playing
E = X(:, 2); N = X(:, 1); Up = -X(:, 3);
span = @(v) max(v) - min(v);
pad = 0.05 * max([span(E), span(N), span(Up), 100]);
axis(scene.AxOverview, [min(E) - pad, max(E) + pad, min(N) - pad, max(N) + pad, ...
    min(0, min(Up)) - pad, max(Up) + pad]);

show(1);
if strcmp(fig.Visible, 'on') && opts.Autoplay
    play(); % returns straight away, the timer does the rest
end
if nargout > 0
    varargout{1} = scene; % (not cleared otherwise, the playback still uses it)
end
end

function [t, X, TV, Um, Ps] = from_timeseries(xs, ts, us, P)
Ps = P;
if ~isstruct(Ps)
    Ps = bsim.rocket.unpack(P);
end
t = xs.Time;
X = rows(xs.Data, 14);
TV = on_times(ts.Time, by_time(ts.Data, numel(ts.Time)), t);
Um = [];
if ~isempty(us)
    % (5 + fins or 9 + fins wide, with the gimbal's motion on the end)
    Um = on_times(us.Time, by_time(us.Data, numel(us.Time)), t);
end
end

% resample onto the state's time steps (holding the last value)
function A = on_times(ta, A, t)
if numel(ta) ~= numel(t) || any(ta ~= t)
    A = interp1(ta, A, t, 'previous', 'extrap');
end
end

function opts = parse(args, u_default)
p = inputParser;
p.addParameter('X', 'xlog');
p.addParameter('Terms', 'terms');
p.addParameter('U', u_default);
p.addParameter('Speed', 1);
p.addParameter('GimbalExaggeration', 1);
p.addParameter('Autoplay', true);
p.parse(args{:});
opts = p.Results;
end

% one row per time step, whatever width the signal is
function A = by_time(D, n_t)
A = squeeze(D);
if size(A, 1) ~= n_t
    A = A.';
end
end

% Timeseries data comes out as N x n or n x 1 x N depending on the signal
function A = rows(D, n)
A = squeeze(D);
if size(A, 1) == n && size(A, 2) ~= n
    A = A.';
end
end
