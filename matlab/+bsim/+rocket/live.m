function live(P, t, x, terms, u, exaggeration)
%LIVE draw the rocket while the simulation runs (call it from a MATLAB Function block)
%   bsim.rocket.live(P, t, x, terms)
%   bsim.rocket.live(P, t, x, terms, u)                    u = the controls, so the fins get drawn
%   bsim.rocket.live(P, t, x, terms, u, gimbal_exaggeration)
%
%   in simulink, a MATLAB Function block with its sample time set on the block (a multiple of the
%   solver step, e.g. set_param(blk, 'SystemSampleTime', '0.05')):
%
%       function show(P, t, x, terms, u)
%       coder.extrinsic('bsim.rocket.live');
%       bsim.rocket.live(P, t, x, terms, u);
%       end
%
%   the window is made on the first call and reused, a new run (t going backwards) starts the
%   trail again. turn on Simulation Pacing (Run > Simulation Pacing) to watch it in real time,
%   otherwise it goes as fast as it can draw. bsim.rocket.replay plays a run back afterwards
persistent scene last_t
if nargin < 5
    u = [];
end
if nargin < 6
    exaggeration = 1;
end

if isempty(scene) || ~isvalid(scene) || ~isvalid(scene.Fig)
    scene = bsim.rocket.Scene(P, 'GimbalExaggeration', exaggeration);
    last_t = -inf;
elseif t < last_t
    scene.reset();
end
scene.GimbalExaggeration = exaggeration;

scene.update(t, x, terms, [], u);
last_t = t;
drawnow limitrate
end
