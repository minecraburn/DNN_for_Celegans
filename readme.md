# code for python

The code is divided into serveral blocks. For use the model, please run (initial module) $\rightarrow$ (user's setting) $\rightarrow$ (loading data module) $\rightarrow$ (loding model module). You can replace the string in (loding model module) to use other module. Note that model in "./mod/re" may result the wrong time label in plot modules.

Then run (curve fit setting module) and then you can run the next four curve fitting module to simulate traces in different conditions
You can change $0$ in
$$ \text{rd=torch.randn(simnu,dim)*(0)} $$
to other values to set the diffusion coeficient, and note that you should uncomment the line
$$\text{x1=x1+rd}.$$

To plot trace, run (the prepare function for ploting) $\rightarrow$ (plot trace in 3D) $\rightarrow$, or 
(the prepare function for ploting 2D) $\rightarrow$ (plot trace in 2D) $\rightarrow$ (plot trace in 2D : color corresponding to behavior). You can change variable "see", "ne1" and "ne2" to plot different traces in different axes. 

Function "learnst(dataofst,tf,use_type=1)" is the function for getting behavior states. use_type=1 means use clustering approach while use_type=2 means use DNN approach.

Module (extra work 1) is for calculating state transition probability in different noise level.
Module (extra work 2) is for calculating state probability for different inhibitions.

For estimate D, please use (esti old D) for model outsides folder "re" and use (esti D) for model insides folder "re".

# code for matlab

1. First load WT_NoStim.mat to load all the neuron activity data, use sorter.m and ImPCA to sort out the origin data, get the PCAs and plot data traces. 

2. To plot common landscape (without inhibition), please do 1 first and then load all file with ".mat" in folder "landscape" and run sigofcy.m. You can change coefficients in "setting field" block to change axes and plotting ares. There are some groups of coefficients in this filed showing how pictures in our article are plotted. If you plot the same figure (same axes) for the second time, you can set the re_la and re_cy1 to 0 to make it quicker. Note to change it back to 1 when you plot other figure (axes).

3. To plot landscape after inhibitions, please do 1 and 2 first and load coe(xx)poi.mat, coe(xx)sig.mat and cycle1coe(xx).mat, cycle2coe(xx).mat (if exists). Then run sigofcydis.m and using the corresponding settings in the "setting field" block.

4. To plot landscape with flux or only flux in the 2-dimesional OLS projection, please do 1 and 2 first (Note that should not do 3, and do 2 with coefficient setting groups "figlandflux"), and load all file with ".mat" in folder "flux". Run fluxland.m / flux.m to plot the figure.

5. For extra1 and extra2, the purpose are same for extra1/2 in python code. Please run 1 and run statetr.m before running e1.m and e2.m. To run e1.m, please first load tranp first and then run e1.m.
To run e2.m, please first load two linep in linep2_1.mat and linep2_2.mat and make linep to be the sum of them.

Update：figures are disordered, thus comment in the code may not match the figure in paper.