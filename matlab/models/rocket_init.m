run(fullfile(fileparts(mfilename('fullpath')), '..', 'setup.m'));
[S, x0, act] = bsim.rocket.params('example');
Prk = bsim.rocket.pack(S);

% open bsim.rocket.replay once a run has finished (the model's StopFcn), false to skip it
show_replay = true;

% --- flight computer settings (the fc struct the Flight Computer blocks read) ---
fc.target = [0; 150];             % where to land [north; east] from the pad [m]
fc.m_cut = S.m_dry + 0.25 * 60;   % MECO with 25% of the propellant left [kg]

% engine and actuator limits
fc.T_vac = S.engine.T_vac_max;    % full throttle vacuum thrust [N]
fc.p_Ae = 101325 * S.engine.A_e;  % thrust lost to sea level pressure [N]
fc.throttle_min = S.engine.throttle_min;
fc.x_G = S.engine.x_G;            % gimbal pivot station, the gimbal's moment arm is x_cm - x_G [m]
fc.gimbal_max = deg2rad(8);
fc.fin_max = deg2rad(20);

% attitude loop (both the gimbal and the fins)
fc.wn = 3;                        % how fast it pulls the nose round [rad/s]
fc.zeta = 0.7;                    % damping
fc.k_roll = 10;                   % how hard it stops roll [1/s]
fc.c_aero = 1.74;                 % aero moment per rad off the airflow per Pa of qbar, from the model

% fin mixer: [roll; pitch; yaw] moment / qbar -> 4 fin angles, measured falling tail first
xd = x0; xd(3) = -1500; xd(4:6) = [0; 0; 100]; xd(14) = fc.m_cut;
uf = [0; 0; 0; Inf; 1; zeros(S.fins.n, 1)]; airf = bsim.rocket.air(1500);
[~, T0] = bsim.rocket.xdot(Prk, xd, uf, airf, 0);
Bfin = zeros(3, S.fins.n);
for k = 1:S.fins.n
    uk = uf; uk(5 + k) = deg2rad(1);
    [~, Tk] = bsim.rocket.xdot(Prk, xd, uk, airf, 0);
    Bfin(:, k) = (Tk.M_body - T0.M_body) / deg2rad(1) / T0.qbar;
end
fc.Mfin = pinv(Bfin);

% steering on the fins (falling tail first)
fc.k_tilt = 0.00211;              % sideways accel per rad of lean per Pa, with the fins holding it
fc.lean_max = deg2rad(3);         % about as far as the fins can hold it off the airflow
fc.kp_fin = 0.05;                 % pull toward the target [1/s^2]
fc.kd_fin = 0.4;                  % damp sideways speed [1/s]

% landing burn
fc.h_ign_max = 2000;              % never relight above this [m]
fc.h_margin = 20;                 % relight early enough to stop this far above the ground [m]
fc.a_ign = 13;                    % relight when stopping needs this much deceleration [m/s^2]
fc.a_land = 12;                   % then slow at this rate (a bit less, to leave margin) [m/s^2]
fc.h_flare = 1;                   % come down the last bit at v_touch [m]
fc.v_touch = 1.2;                 % touchdown speed [m/s]
fc.kv = 2.5;                      % vertical speed loop [1/s]
fc.tilt_max = 0.12;               % sideways accel up to 0.12 of the vertical (about 7 deg)
fc.t_final = 2;                   % last part of the burn: just stop sliding and stand up [s]
fc.tau_final = 0.8;               % how quickly it kills sideways speed then [s]
fc.h_cutoff = 0.15;               % engine off below this [m]
