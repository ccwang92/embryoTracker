# GPU-accelerated single-pass raycaster

This project is a simple visualiser based on GPU-accelerated single-pass
[volumetric raycasting](https://en.wikipedia.org/wiki/Volume_ray_casting),
implemented in
[GLSL](https://www.khronos.org/opengl/wiki/OpenGL_Shading_Language) and C++
with the [Qt](https://www.qt.io/) framework. It aims to provide a basic
skeleton for a visualiser that can be easily extended with a feature-rich GUI.
A brief overview of the design and implementation is described in a [blog
post](https://martinopilia.com/posts/2018/09/17/volume-raycasting.html).

Three simple examples of shaders are provided for
[isosurface](https://en.wikipedia.org/wiki/Isosurface) rendering, front-to-back
[alpha-blending](https://en.wikipedia.org/wiki/Alpha_compositing), and [maximum
intensity
projection](https://en.wikipedia.org/wiki/Maximum_intensity_projection). It
implements a simple reader for the [VTK Structured Point
legacy](http://www.cs.utah.edu/~ssingla/Research/file-formats.pdf) file format,
allowing to load volumes from file out-of-the-box.


# License

The software is distributed under the MIT license.
