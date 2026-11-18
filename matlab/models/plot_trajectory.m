out = sim('rocket_plant');
set_param(gcb, 'SystemSampleTime', '1/30')
bsim.rocket.replay(out, Prk)                          
bsim.rocket.replay(out, Prk, 'GimbalExaggeration', 3)