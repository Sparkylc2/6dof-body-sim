# 6dof-body-sim
A 6 dof sim with the core based on Small Unmanned Aircraft: Theory and Practice by Beard and McLain.

The notes attached show my own working while looking at the text book, and functions implementing stuff from the textbook will have relevant sections reference in the comments (for things that come from the textbook).

Currently, the body model, wind models, sensor models, and basic control piping have been set up, alongside a satellite model that currently includes the J2 equator bulge, gravity gradients, and reaction wheels (alongside thrusters, and variable COM based on fuel loss).


There is now also a small polynomial SINDy model, based on the tutorial repo and the concepts explored in this [paper](https://arxiv.org/pdf/2607.15077).

I hope to have the time to expand this further and follow more parts of the textbook + my ideas, continuing into full autopilot design and state estimation.
