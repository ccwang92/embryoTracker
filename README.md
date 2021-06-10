# A unified cell segmentation, tracking and visualization platform
This is an on-going project aiming to build a platform that supports 3D+time (colorful or gray) imaging data visualization, cell segmentation, cell tracking and results curation. It is built using Qt with C++ and OpenGL.

## Current interface and visualization quality
<p align="center">
  <img height="500" src="imgs/org_im3d.png">
</p>

## Current segmentation quality
<p align="center">
  <img height="500" src="imgs/seg_res.png">
</p>


# Run the software
The software is still under development. The code was tested with Ubuntu 18.04.5 LTS and Ubuntu 16.04 LTS. The IDE and compiler are Qt 5.15.2 and GCC 5.3.1. Newer Linux, Qt, or GCC versions should also work.
1. Clone the project to your local hard drive.

2. Download the 3rd-party vaa3d codes [here](https://drive.google.com/drive/folders/18ZK8gdnF8-MK2oPed-6jtIGr6tMEXAia?usp=sharing) for tiff data loading and put it (the folder 'src_3rd') at the root folder of the project.

3. Install OpenCV and Boost libraries to your system.

4. Use Qt to compile the code by opening the .pro file.


# Other resources
[synQuant](https://github.com/yu-lab-vt/SynQuant), [muSSP](https://github.com/yu-lab-vt/muSSP), [CINDA](https://github.com/yu-lab-vt/CINDA), [Vaa3d project.](https://github.com/Vaa3D), [RayCasting volume rendering](https://martinopilia.com/posts/2018/09/17/volume-raycasting.html) and [QDarkStyleSheet.](https://github.com/ColinDuquesnoy/QDarkStyleSheet)

# Reference

[1] Congchao Wang, Yizhi Wang, Guoqiang Yu, Efficient Global Multi-object Tracking Under Minimum-cost Circulation Framework, arXiv:1911.00796. (accepted by IEEE Trans. on PAMI)

[2] Congchao Wang, Yizhi Wang, Yinxue Wang, Chiung-Ting Wu, Guoqiang Yu. (2019). muSSP: Efficient min-cost flow algorithm for multi-object tracking. In Advances in Neural Information Processing Systems (pp. 425-434).

[3] Yizhi Wang*, Congchao Wang*, Petter Ranefall, Gerard Joey Broussard, Yinxue Wang, Guilai Shi, Boyu Lyu, Chiung-Ting Wu, Yue Wang, Lin Tian, Guoqiang Yu. (2020). SynQuant: An Automatic Tool to Quantify Synapses from Microscopy Images, Bioinformatics, 36(5), 1599–1606.

[4] Peng, H., Ruan, Z., Long, F., Simpson, J.H., and Myers, E.W. (2010) “V3D enables real-time 3D visualization and quantitative analysis of large-scale biological image data sets,” Nature Biotechnology, Vol. 28, No. 4, pp. 348-353.

[5] Peng, H, Ruan, Z., Atasoy, D., and Sternson, S. (2010) “Automatic reconstruction of 3D neuron structures using a graph-augmented deformable model,” Bioinformatics, Vol. 26, pp. i38-i46, 2010.
