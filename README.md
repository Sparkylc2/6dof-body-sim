# 6dof-body-sim
A 6 dof sim with the core based on Small Unmanned Aircraft: Theory and Practice by Beard and McLain.

The notes attached show my own working while looking at the text book, and functions implementing stuff from the textbook will have relevant sections reference in the comments (for things that come from the textbook).

### Notes index (`six dof sim.pdf`)
- Frames & rotations (Ch 2) - p2
- Kinematics & rigid-body dynamics (Ch 3) - p3
- Quaternion dynamics - p4
- Forces & moments, wind/turbulence (Ch 4) - p5
- Sensor models (Ch 7) - p6 – p7
- Basic PID control - p7
- Spacecraft model (orbit, J2, gravity gradient, reaction wheels, thrusters) - p8
- State estimation (observers, Kalman filter / EKF) - p9
- Spacecraft attitude control (reaction wheel pointing) - p10
- SINDy - p11 – p12
- Successive convexification - p13 - p14

Currently, the body model, wind models, sensor models, and basic control piping have been set up, alongside a satellite model that currently includes the J2 equator bulge, gravity gradients, and reaction wheels (alongside thrusters, and variable COM based on fuel loss).


There is now also a small polynomial SINDy model, based on the tutorial repo and the concepts explored in this [paper](https://arxiv.org/pdf/2607.15077). This was used to create a small fit on the phugoid mode of an aircraft (this is a particularly simple example I was just curious as to how the method worked).

In lieu of work I'm doing with Imperials rocketry society, implementing a successive convexification scheme (specifically for the self landing goals of the society) will be included! The initial implementation will be based on this [paper](https://arxiv.org/abs/2405.00061). Regarding the core conic/MPC solvers I will probably have them AI generated for now (I don't feel like adding a solver dependency) 

I hope to have the time to expand this further and follow more parts of the textbook + my ideas, continuing into full autopilot design and state estimation.
