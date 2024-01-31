Download entire BreastDCE_DRO repo and run 'example.m'

example.m
-----------------------------------------------------------
1. Run Data collection to read all data

2. Generate DRO:

   2-1. Without any option, a randomly selected case (benign or malignant) will be generated with randomly selected PK parameters and AIFs.

   2-2. For external input of parMap: user can supply their own parameters.
   The parameter maps need to be in 2D, and need to be generated using the segmentation masks.
   Glandular, muscle and skin need to have 3 parameters [0, vp, Fp, PS], while rest regions need to use 4 parameters for TCM (i.e. [ve, vp, Fp, PS])

   2-3. For external input of B1 maps: need to supply as an option (i.e. option.B1 = B1map)

   2-4. For external input of AIF: need to supply as an option (i.e. option.AIF = aif)

3. Radial k-space data generation: Current implementaion generates golden-angle separated radial trajectories based on the desired number of spokes per each temporal frame for sampling.

4. GRASP recon using BART: Need BART toolkit installed for running recon. Please include BART directory in the path. Current command runs BART recon with TV regularization (lamb=0.01).

   
