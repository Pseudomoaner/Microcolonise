# Microcolonise

Analysis scripts associated with the paper 'Droplet printing reveals the importance of micron-scale structure for bacterial ecology'. The scripts are split into two folders:

## 2DAnalyses

Analyses that allow the automated registration of prints to the x-y coordinate space, removal of microcolonies far from the main print, generation of radial profiles and measurement of 2D microcolony features. To run, use the script `Do2DAnalyses.m`.

Please note that these analyses rely on performing segmentation and feature extraction with [FAST](https://github.com/Pseudomoaner/FAST).

## 3DAnalyses

Analyses that allow 3D segmentation of microcolony data. For full details of the segmentation approach used, please refer to the publications in the references section. To run, use the script `run3Dsegmentation.m`.

## References

- Jin, J., Yang, L., Zhang, X. & Ding, M. (2013). Vascular tree segmentation in medical images using Hessian-based multiscale filtering and level set method. Computational and Mathematical Methods in Medicine, 502013. https://doi.org/10.1155/2013/502013
- Kumar, R.K., Meiller-Legrand, T., Alcinesio, A., Gonzalez, D., Mavridou, D.A.I., Meacock, O.J., Smith, W.P.J., Zhou, L., Kim, W., Pulcu, G., Bayley, H. & Foster, K.R. (2021). Droplet printing reveals the importance of micron-scale structure for bacterial ecology. Nat Commun **12**, 857. https://doi.org/10.1038/s41467-021-20996-w
