/*
 * Copyright © 2018 Martino Pilia <martino.pilia@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


#include "raycastvolume.h"
//#include "vtkvolume.h"
#include <QRegularExpression>
#include "raycastcanvas.h"
#include <algorithm>
#include <cmath>


/*!
 * \brief Create a two-unit cube mesh as the bounding box for the volume.
 */
RayCastVolume::RayCastVolume(void* _glwidget)
    : m_volume_texture {0}
    , m_noise_texture {0}
    , m_cube_vao {
          {
              -1.0f, -1.0f,  1.0f,
               1.0f, -1.0f,  1.0f,
               1.0f,  1.0f,  1.0f,
              -1.0f,  1.0f,  1.0f,
              -1.0f, -1.0f, -1.0f,
               1.0f, -1.0f, -1.0f,
               1.0f,  1.0f, -1.0f,
              -1.0f,  1.0f, -1.0f,
          },
          {
              // front
              0, 1, 2,
              0, 2, 3,
              // right
              1, 5, 6,
              1, 6, 2,
              // back
              5, 4, 7,
              5, 7, 6,
              // left
              4, 0, 3,
              4, 3, 7,
              // top
              2, 6, 7,
              2, 7, 3,
              // bottom
              4, 5, 1,
              4, 1, 0,
          }
      }
{
    initializeOpenGLFunctions();
    glWidget = _glwidget;
}


/*!
 * \brief Destructor.
 */
RayCastVolume::~RayCastVolume()
{
}


/*!
 * \brief get the volume from main window.
 * \param frame4show starts from 0.
 */
void RayCastVolume::transfer_volume(void*data, double p_min, double p_max, long sx,
                                    long sy, long sz, long sc) {

    if (sc < 1) //gray image that stacked by channel, sc will equal to time
    {
        throw std::runtime_error("data to show is NULL.");
    }
//    VTKVolume volume {filename.toStdString()};
//    volume.uint8_normalised();
//    m_size = QVector3D(std::get<0>(volume.size()), std::get<1>(volume.size()), std::get<2>(volume.size()));
//    m_origin = QVector3D(std::get<0>(volume.origin()), std::get<1>(volume.origin()), std::get<2>(volume.origin()));
//    m_spacing = QVector3D(std::get<0>(volume.spacing()), std::get<1>(volume.spacing()), std::get<2>(volume.spacing()));
//    m_range = volume.range();
//    data = volume.data();

    //this->set_size(QVector3D(sx, sy, sz));
    //this->set_origin(QVector3D(0,0,0));
    //this->set_spacing(QVector3D(1,1,1)); // default all one per voxel
    //this->set_range(p_min, p_max);
    m_size = QVector3D(sx, sy, sz);
    m_origin = QVector3D(0,0,0);
    m_spacing = QVector3D(1,1,1);
    m_range = std::make_pair(p_min, p_max);

    glDeleteTextures(1, &m_volume_texture);
    glGenTextures(1, &m_volume_texture);
    glBindTexture(GL_TEXTURE_3D, m_volume_texture);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  // The array on the host has 1 byte alignment

    glTexImage3D(GL_TEXTURE_3D, 0, GL_R8, m_size.x(), m_size.y(), m_size.z(),
                 0, GL_RED, GL_UNSIGNED_BYTE, data);
    glBindTexture(GL_TEXTURE_3D, 0);

    // initial values for bounding box
    boundingBox.x0 = boundingBox.y0 = boundingBox.z0 = 0.0;
    boundingBox.x1 = sx * m_spacing[0];
    boundingBox.y1 = sy * m_spacing[1];
    boundingBox.z1 = sz * m_spacing[2];
}


/*!
 * \brief Create a noise texture with the size of the viewport.
 */
void RayCastVolume::create_noise(void)
{
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    const int width = viewport[2];
    const int height = viewport[3];

    std::srand(std::time(NULL));
    unsigned char noise[width * height];

    for (unsigned char *p = noise; p <= noise + width * height; ++p) {
        *p = std::rand() % 256;
    }

    glDeleteTextures(1, &m_noise_texture);
    glGenTextures(1, &m_noise_texture);
    glBindTexture(GL_TEXTURE_2D, m_noise_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, width, height, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, noise);
    glBindTexture(GL_TEXTURE_2D, 0);
}


/*!
 * \brief Render the bounding box.
 */
void RayCastVolume::paint(void)
{
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_3D, m_volume_texture);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, m_noise_texture);

    m_cube_vao.paint();

    if (!m_size.isNull())//(!(((RayCastCanvas*)glWidget)->getDataImporter()))
    {
        //((RayCastCanvas*)glWidget)->renderText(0,0,0, "X");
//        if (bShowBoundingBox || bShowAxes || bShowXYTranslateArrows)// a frame box [-1,+1]^3
//        {
//            setBoundingBoxSpace(boundingBox);
//            drawBoundingBoxAndAxes(boundingBox);
//        }
    }
}


/*!
 * \brief Range of the image, in intensity value.
 * \return A pair, holding <minimum, maximum>.
 */
std::pair<double, double> RayCastVolume::range() {
    return m_range;
}


/*!
 * \brief Scale factor to model space.
 *
 * Scale the bounding box such that the longest side equals 1.
 */
float RayCastVolume::scale_factor(void)
{
    auto e = m_size * m_spacing;
    return std::max({e.x(), e.y(), e.z()});
}


/*** The following functions are for bounding box and axes draw ****/
void RayCastVolume::setBoundingBoxSpace(BoundingBox BB)
{
    float DX = BB.Dx();
    float DY = BB.Dy();
    float DZ = BB.Dz();
    float maxD = BB.Dmax();

    double s[3];
    s[0] = 1/maxD *2;
    s[1] = 1/maxD *2;
    s[2] = 1/maxD *2;
    double t[3];
    t[0] = -BB.x0 -DX /2;
    t[1] = -BB.y0 -DY /2;
    t[2] = -BB.z0 -DZ /2;

    // from boundingBox space ==> fit in [-1, +1]^3
    glScaled(s[0], s[1], s[2]);
    glTranslated(t[0], t[1], t[2]);
}

inline void box_quads(const BoundingBox & BB)
{
#define BB_VERTEX(xi,yi,zi)  glVertex3d(BB.x##xi, BB.y##yi, BB.z##zi)

    BB_VERTEX(0, 0, 0);	BB_VERTEX(0, 1, 0);	BB_VERTEX(1, 1, 0); BB_VERTEX(1, 0, 0); //z=0
    BB_VERTEX(0, 0, 1);	BB_VERTEX(0, 1, 1);	BB_VERTEX(1, 1, 1); BB_VERTEX(1, 0, 1); //z=1

    BB_VERTEX(0, 0, 0);	BB_VERTEX(1, 0, 0);	BB_VERTEX(1, 0, 1); BB_VERTEX(0, 0, 1); //y=0
    BB_VERTEX(0, 1, 0);	BB_VERTEX(1, 1, 0);	BB_VERTEX(1, 1, 1); BB_VERTEX(0, 1, 1); //y=1

    BB_VERTEX(0, 0, 0);	BB_VERTEX(0, 0, 1);	BB_VERTEX(0, 1, 1); BB_VERTEX(0, 1, 0); //x=0
    BB_VERTEX(1, 0, 0);	BB_VERTEX(1, 0, 1);	BB_VERTEX(1, 1, 1); BB_VERTEX(1, 1, 0); //x=1

}

inline void draw_tri(const XYZ P1, const XYZ P2, const XYZ P3, const XYZ offst)
{
    // front/back
    glBegin(GL_TRIANGLES);
    glVertex3f(P1.x, P1.y, P1.z);
    glVertex3f(P2.x, P2.y, P2.z);
    glVertex3f(P3.x, P3.y, P3.z);
    glVertex3f(P1.x+offst.x, P1.y+offst.y, P1.z+offst.z);
    glVertex3f(P2.x+offst.x, P2.y+offst.y, P2.z+offst.z);
    glVertex3f(P3.x+offst.x, P3.y+offst.y, P3.z+offst.z);
    glEnd();
    // sides
    glBegin(GL_QUADS);
    glVertex3f(P1.x, P1.y, P1.z);
    glVertex3f(P1.x+offst.x, P1.y+offst.y, P1.z+offst.z);
    glVertex3f(P2.x+offst.x, P2.y+offst.y, P2.z+offst.z);
    glVertex3f(P2.x, P2.y, P2.z);
    glVertex3f(P2.x, P2.y, P2.z);
    glVertex3f(P2.x+offst.x, P2.y+offst.y, P2.z+offst.z);
    glVertex3f(P3.x+offst.x, P3.y+offst.y, P3.z+offst.z);
    glVertex3f(P3.x, P3.y, P3.z);
    glVertex3f(P3.x, P3.y, P3.z);
    glVertex3f(P3.x+offst.x, P3.y+offst.y, P3.z+offst.z);
    glVertex3f(P1.x+offst.x, P1.y+offst.y, P1.z+offst.z);
    glVertex3f(P1.x, P1.y, P1.z);
    glEnd();
}

void RayCastVolume::drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth, float AlineWidth)
{
    glPushAttrib(GL_LINE_BIT | GL_POLYGON_BIT);
            //| GL_DEPTH_BUFFER_BIT);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_POLYGON_OFFSET_LINE);

//	glPolygonOffset(0, -1); // deal z-fighting, 081120
//	glDepthFunc(GL_LEQUAL);
    if (posXTranslateBB != 0) delete posXTranslateBB;
    if (negXTranslateBB != 0) delete negXTranslateBB;
    if (posYTranslateBB != 0) delete posYTranslateBB;
    if (negYTranslateBB != 0) delete negYTranslateBB;
    posXTranslateBB=0, negXTranslateBB=0,
            posYTranslateBB=0, negYTranslateBB=0;
    // an indicator of coordinate direction
    if (bShowAxes && AlineWidth>0)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        float D = (BB.Dmax());
        float ld = D*0.0001; //1e-4 is best
        float td = D*0.015;
        XYZ A0 = BB.Vabsmin();
        XYZ A1 = BB.V1() + D*0.05;

        glPolygonOffset(-0.002, -2); //(-0.002, -2) for good z-fighting with bounding box, 081120,100823

        glLineWidth(AlineWidth); // work only before glBegin(), by RZC 080827
        //glBegin(GL_QUADS);
//        glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
//        {
////            glColor3f(1, 0, 0);		box_quads( BoundingBox(A0, XYZ(A1.x, A0.y+ld, A0.z+ld)) );
////            glColor3f(0, 1, 0);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A1.y, A0.z+ld)) );
////            glColor3f(0, 0, 1);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A0.y+ld, A1.z)) );
//            glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(10, 0, 0);
//            glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 10, 0);
//            glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 10);
//        }
//        glEnd();

        glColor3f(1, 0, 0);		drawString(A1.x+td, A0.y, A0.z, "X", 1, 0);
        glColor3f(0, 1, 0);		drawString(A0.x, A1.y+td, A0.z, "Y", 1, 0);
        glColor3f(0, 0, 1);		drawString(A0.x, A0.y, A1.z+td, "Z", 1, 0);
    }

    if (bShowBoundingBox && BlineWidth>0)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        glPolygonOffset(0, -1); // deal z-fighting with volume, 081120

        glLineWidth(BlineWidth); // work only before glBegin(), by RZC 080827
        glBegin(GL_QUADS);
        //glBegin(GL_LINES);
        {
            glColor3fv(color_line.c);	box_quads(BB);
        }
        glEnd();
    }

    glPopAttrib();
}
void RayCastVolume::drawString(float x, float y, float z, const char* text, int shadow, int fontsize)
{
    if (! glWidget)  return;

    if (shadow)
    {
        glPushAttrib(GL_DEPTH_BUFFER_BIT);
        glPushAttrib(GL_CURRENT_BIT);
            //glColor3ub(50,50,50);
            //glColor3ub(200,200,200);

            // CMB MSVC debugger with Qt 4.7 triggers assert if font weight > 99
            // QFont f;  f.setPointSize(f.pointSize()+1); f.setWeight(f.weight()+200);
            QFont f;  f.setPointSize(f.pointSize()+1); f.setWeight(99);
//#if defined(USE_Qt5)
//#else
//            ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text), f);
//#endif
            QPainter painter((QOpenGLWidget*)glWidget);
            painter.setPen(QColor(50,50,50));
            painter.setFont(QFont("Arial", 16));
            painter.drawText(QPoint(int(200),int(200)), QString(text));
            painter.end();
        glPopAttrib();
        glDepthFunc(GL_LEQUAL);
    }

    QFont f1;  f1.setPointSize((fontsize>0)?fontsize:30); f1.setWeight(99);
//    if (fontsize>0)
//#if defined(USE_Qt5)
//#else
//        ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text), f1);
//    else
//        ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text));
//#endif


    if (shadow)
    {
        glPopAttrib();
    }
}
