This repository holds code used to generate the figures and tables in our paper "Sensitivity to control signals in triphasic rhythmic neural systems: a comparative mechanistic analysis via infinitesimal local timing response curves". Figures 1-4, 7-9, 11, 12, and 14 are for relaxation-oscillator model. Figures 5, 15, and 16 are for the heteroclinic cycling model. Figures 6, 17, and 18 are for the competitive threshold-linear model. The rest figures are not generated via simulation.

Fig. 1 shows a typical solution of intrinsic-release system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_1.m

Fig. 2 shows a typical solution of intrinsic-escape system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.24, 'E', 0.14); Theta = struct('h', -39, 'mp', -37); Sigma = struct('h', 9, 'mp', -6);", and then run Fig_2.m

Fig. 3 shows a typical solution of synaptic-release system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_3.m

Fig. 4 shows a typical solution of synaptic-escape system. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 5, 'mp', -6);" in relaxation.m, and then run Fig_4.m

Fig. 5 shows a heteroclinic cycling solution and a limit cycle solutions of piecewise linear Aplysia model. To generate figure, run Fig_5.m

Fig. 6 shows a solution of competitive threshold-linear model. To generate figure, run Fig_6.m

Fig. 7 shows phase plane for cell 1 in intrinsic release, with several values of d1. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_7.m

Fig. 8 shows phase plane for cell 1 in synaptic release, with several values of d1. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_8.m

Fig. 9 shows lTRC for (a) intrinsic release and (b) synaptic release. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_9a and Fig_9b.

Fig. 11 shows lTRC for synaptic escape. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 5, 'mp', -6);" in relaxation.m, and change the event functions to be "value=P(3)-theta_imag". Then run Fig_11.m

Fig. 12 shows phase plane for three cells in synaptic escape, with several values of d1. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 5, 'mp', -6);" in relaxation.m, and then run Fig_12.m

Fig. 14 shows a double-period solution of intrinsic-escape system.  To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.24, 'E', 0.14); Theta = struct('h', -39, 'mp', -37); Sigma = struct('h', 9, 'mp', -6);", and then run Fig_14.m

Fig. 15 shows solutions of piecewise linear Aplysia model with several values of a1. To generate figure, run Fig_15.m

Fig. 16 shows lTRC for heteroclinic cycling model. To generate figure, run Fig_16.m

Fig. 17 shows solutions of competitive threshold-linear model in the perturbed case and unperturbed case. To generate figure, run Fig_17.m

Fig. 18 shows lTRC for competitive threshold-linear model. To generate figure, run Fig_18.m

Fig. 19 compares lTRC and iPRC of (a) intrinsic release and (b) synaptic release. To generate figure, change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m, and then run Fig_19.m

Table 1 lists duration changes in intrinsic release. Change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m. To generate data in top line, run Fig_9a.m. To generate data in bottom line, change "delta = -0.05; initialsp = [-10.0000  -62.4740  -63.8739    0.4056    0.7728    0.3977];" in Fig_9a.m and run it

Table 2 lists duration changes in synaptic release. Change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 6, 'mp', -6);" in relaxation.m. To generate data in top line, run Fig_9b.m. To generate data in bottom line, change "delta = -0.05; initialsp = [-10.0000  -62.7965  -63.8955    0.4056    0.7028    0.3903];" in Fig_9b.m and run it

Table 3 lists duration changes in synaptic escape. Change the parameter values "G = struct('NaP', 6.8, 'L', 3, 'I', 0.4, 'E', 0.1); Theta = struct('h', -40, 'mp', -37); Sigma = struct('h', 5, 'mp', -6);" in relaxation.m. Change the event functions to be "value=P(3)-theta_imag". To generate data in top line, run Fig_11.m. To generate data in bottom line, change "delta = -0.01; initialsp = [-10.0000  -62.5848  -63.9052    0.4049    0.7501    0.3878];" in Fig_11.m and run it

Table 4 lists duration changes in heteroclinic cycling model. To generate data in top line, run Fig_16.m. To generate data in bottom line, change "ap = -0.0095; initialsp = [0.500329161858328 0.0162116170731029 0.239209693626476];" in Fig_16.m and run it

Table 5 lists duration changes in competitive threshold-linear model. To generate data in top line, run Fig_18.m. To generate data in bottom line, change "thetap = 0.99; initialsp = [0.500057726382846,0.0122982618497747,0.425056979349457];" in Fig_18.m and run it
