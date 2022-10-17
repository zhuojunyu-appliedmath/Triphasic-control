This repository holds code used to generate the figures and tables in our paper "Control in triphasic rhythmic neural systems: a comparative mechanistic analysis via infinitesimal local timing response curves". Figures 1-4, 7-9, 11, 12, and 14 are for relaxation-oscillator model. Figures 5, 15, and 16 are for the heteroclinic cycling model. Figures 6, 17, and 18 are for the competitive threshold-linear model. The rest figures are not generated via simulation.

Fig. 1 shows a typical solution of intrinsic-release system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_1.m

Fig. 2 shows a typical solution of intrinsic-escape system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.24, 'E', 0.14); Theta = struct('h', -39, 'mp', -37); Sigma = struct('h', 9, 'mp', -6);", and then run Fig_2.m

Fig. 3 shows a typical solution of synaptic-release system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_3.m

Fig. 4 shows a typical solution of synaptic-escape system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 5, 'mp', -6);" in relaxation.m, and then run Fig_4.m

Fig. 5 shows a heteroclinic cycling solution and a limit cycle solutions of piecewise linear Aplysia model. To generate figure, run Fig_5.m

Fig. 6 shows a solution of competitive threshold-linear model. To generate figure, run Fig_6.m

Fig. 7 shows phase plane for cell 1 in intrinsic release, with several values of d1. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_7.m

Fig. 8 shows phase plane for cell 1 in synaptic release, with several values of d1. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_8.m
