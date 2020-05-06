A HYBRID GLOBAL IMAGE ORIENTATION METHOD FOR SIMULTANEOUSLY ESTIMATING GLOBAL ROTATIONS AND GLOBAL TRANSLATIONS

Xin Wang*, Teng Xiao, Yoni Kasten

## Properties:

The code reconstructs images' exterior poses by using pairwise essential matrices and tie points.
our code is only licensed for noncommerical research use, you can cite:


To be appear In Proceedings{
	author = Xin, Wang and Teng Xiao, Yoni Kasten},
	title = {A HYBRID GLOBAL IMAGE ORIENTATION METHOD FOR SIMULTANEOUSLY ESTIMATING GLOBAL ROTATIONS AND GLOBAL TRANSLATIONS},
	booktitle = {XXIV ISPRS Congress 2020},
	month = {June},
	year = {2020}
}

## Installation:

    This work is well tested on the matlab2016b, higher version of Matlab also works

## Use:

   We only computer exterior parameters, you can run run_building.m file on your matlab package, and the result will be visually shown and saved as a txt file.
   Note that, the input are generate by ProcessScript.m, one can generated the necessaries from OPENMVG, which is also used for refining our result.
   The ground truth exterior parameters of fountaion Herzjezu castle are provided by the corresponding publisher 

## Acknowledgement:

We develop our code based on https://github.com/amnonge/GPSFM-code, you may also interesed in those works:

InProceedings{Kasten_2019_ICCV,
	author = {Kasten, Yoni and Geifman, Amnon and Galun, Meirav and Basri, Ronen},
	title = {Algebraic Characterization of Essential Matrices and Their Averaging in Multiview Settings},
	booktitle = {The IEEE International Conference on Computer Vision (ICCV)},
	month = {October},
	year = {2019}
}

InProceedings{Kasten_2019_CVPR,
author = {Kasten, Yoni and Geifman, Amnon and Galun, Meirav and Basri, Ronen},
title = {GPSfM: Global Projective SFM Using Algebraic Constraints on Multi-View Fundamental Matrices},
booktitle = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
month = {June},
year = {2019}
}
