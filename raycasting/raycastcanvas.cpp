#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>

#include <QtWidgets>

#include "raycastcanvas.h"


/*!
 * \brief Convert a QColor to a QVector3D.
 * \return A QVector3D holding a RGB representation of the colour.
 */
QVector3D to_vector3d(const QColor& colour) {
    return QVector3D(colour.redF(), colour.greenF(), colour.blueF());
}


/*!
 * \brief Constructor for the canvas.
 * \param parent Parent widget.
 */
RayCastCanvas::RayCastCanvas(QWidget *parent)
    : QOpenGLWidget {parent},
    m_stepLength(0.003),
    m_background(QColor(204, 255, 255)),
    m_raycasting_volume {nullptr},
    m_active_mode("MIP")
{
    // Register the rendering modes here, so they are available to the UI when it is initialised
    m_modes["Isosurface"] = [&]() { RayCastCanvas::raycasting("Isosurface"); };
    m_modes["Alpha blending"] = [&]() { RayCastCanvas::raycasting("Alpha blending"); };
    m_modes["MIP"] = [&]() { RayCastCanvas::raycasting("MIP"); };

    // set focus policy to accept key press events
    this->setFocusPolicy(Qt::StrongFocus);
}


/*!
 * \brief Destructor.
 */
RayCastCanvas::~RayCastCanvas()
{
    for (auto& [key, val] : m_shaders) {
        delete val;
    }
    delete m_raycasting_volume;
}


/*!
 * \brief Initialise OpenGL-related state.
 */
void RayCastCanvas::initializeGL()
{
    initializeOpenGLFunctions();

    if(hasMouseTracking())
        setMouseTracking(false);
    //canvas_painter = new QPainter(this);
    m_raycasting_volume = new RayCastVolume(this);
    //m_raycasting_volume->create_noise();

    add_shader("Isosurface", ":/shaders/isosurface.vert", ":/shaders/isosurface.frag");
    add_shader("Alpha blending", ":/shaders/alpha_blending.vert", ":/shaders/alpha_blending.frag");
    add_shader("MIP", ":/shaders/maximum_intensity_projection.vert", ":/shaders/maximum_intensity_projection.frag");
    //glClearColor(1.f, 1.f, 1.f, 1.f);
}


/*!
 * \brief Callback to handle canvas resizing.
 * \param w New width.
 * \param h New height.
 */
void RayCastCanvas::resizeGL(int w, int h)
{
    (void) w; (void) h;
    m_viewportSize = {(float) scaled_width(), (float) scaled_height()};
    m_aspectRatio = (float) scaled_width() / scaled_height();
    glViewport(0, 0, scaled_width(), scaled_height());
    //glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    //m_raycasting_volume->create_noise();
}


/*!
 * \brief Paint a frame on the canvas.
 */
void RayCastCanvas::paintGL()
{
//    glClear(GL_COLOR_BUFFER_BIT);

//    QPainter canvas_painter(this);
//    canvas_painter.setPen(Qt::black);
//    canvas_painter.setFont(QFont("Arial", 16));
//    canvas_painter.drawText(0, 0, scaled_width(), scaled_height(), Qt::AlignCenter, "Hello World!");
//    canvas_painter.end();

    // Compute geometry
    m_viewMatrix.setToIdentity();
    //m_viewMatrix.translate(centerShift->x()/width(), centerShift->y()/height(), -4.0f * std::exp(m_distExp / 600.0f));
    m_viewMatrix.translate(centerShift->x(), centerShift->y(), -4.0f * std::exp(m_distExp / 600.0f));
    m_viewMatrix.rotate(m_trackBall.rotation());

    m_modelViewProjectionMatrix.setToIdentity();
    m_modelViewProjectionMatrix.perspective(m_fov, (float)scaled_width()/scaled_height(), 0.1f, 100.0f);
    m_modelViewProjectionMatrix *= m_viewMatrix * m_raycasting_volume->modelMatrix();

    m_normalMatrix = (m_viewMatrix * m_raycasting_volume->modelMatrix()).normalMatrix();
    //m_normalMatrix = (m_raycasting_volume->modelMatrix()).normalMatrix();
    m_rayOrigin = m_viewMatrix.inverted() * QVector3D({0.0, 0.0, 0.0});
    //m_rayOrigin = QVector3D({0.0, 0.0, 0.0});
    // Perform raycasting
    m_modes[m_active_mode]();

//    qDebug("The canvas size is : %d, %d, and the volume size is : %f, %f %f\n",
//           this->scaled_width(), this->scaled_height(), //this->scaled_depth(),
//           this->m_raycasting_volume->get_size().x(), this->m_raycasting_volume->get_size().y(),
//           this->m_raycasting_volume->get_size().z());

    if (bShowBoundingBox)
    {
        glPushMatrix();
        draw_bbox();
        glPopMatrix();
    }
}


/*!
 * \brief Width scaled by the pixel ratio (for HiDPI devices).
 */
GLuint RayCastCanvas::scaled_width()
{
    return devicePixelRatio() * width();
}


/*!
 * \brief Height scaled by the pixel ratio (for HiDPI devices).
 */
GLuint RayCastCanvas::scaled_height()
{
    return devicePixelRatio() * height();
}


/*!
 * \brief Perform isosurface raycasting.
 */
void RayCastCanvas::raycasting(const QString& shader)
{
    m_shaders[shader]->bind();
    {
        m_shaders[shader]->setUniformValue("ViewMatrix", m_viewMatrix);
        m_shaders[shader]->setUniformValue("ModelViewProjectionMatrix", m_modelViewProjectionMatrix);
        m_shaders[shader]->setUniformValue("NormalMatrix", m_normalMatrix);
        m_shaders[shader]->setUniformValue("aspect_ratio", m_aspectRatio);
        m_shaders[shader]->setUniformValue("focal_length", m_focalLength);
        m_shaders[shader]->setUniformValue("viewport_size", m_viewportSize);
        m_shaders[shader]->setUniformValue("ray_origin", m_rayOrigin);
        m_shaders[shader]->setUniformValue("top", m_raycasting_volume->top());
        m_shaders[shader]->setUniformValue("bottom", m_raycasting_volume->bottom());
        m_shaders[shader]->setUniformValue("background_colour", to_vector3d(m_background));
        m_shaders[shader]->setUniformValue("light_position", m_lightPosition);
        m_shaders[shader]->setUniformValue("material_colour", m_diffuseMaterial);
        m_shaders[shader]->setUniformValue("step_length", m_stepLength);
        m_shaders[shader]->setUniformValue("threshold", m_threshold);
        m_shaders[shader]->setUniformValue("gamma", m_gamma);
        m_shaders[shader]->setUniformValue("volume", 0);
        m_shaders[shader]->setUniformValue("jitter", 1);
        m_shaders[shader]->setUniformValue("consider_transparency", consider_transparency /*0*/);

        glClearColor(m_background.redF(), m_background.greenF(), m_background.blueF(), m_background.alphaF());
        glClear(GL_COLOR_BUFFER_BIT);

        m_raycasting_volume->paint();
    }
    m_shaders[shader]->release();
}
/*!
 * \brief import the data and set it as color RGBA data if needed
 */
void RayCastCanvas::setVolume(long frame4display) {
    if (!data_importer)
    {
        throw std::runtime_error("data_importer has not been initialized.");
    }
    if (!data_importer->p_vmin){// if max min value not defined
        data_importer->updateminmaxvalues();
    }
    if (!total_rgbaBuf) // if buffer has not been initialized
    {
        QElapsedTimer timer;
        try
        {
            if (data_importer->image4d && data_importer->image4d->getCDim()>0)
            {
                //              data_unitbytes = data_importer->image4d->getUnitBytes();
                //              data4dp = data_importer->image4d->getRawData();
                //              data4d_uint8 = data_importer->image4d->data4d_uint8;

                size1= dim1 = data_importer->image4d->getXDim();
                size2= dim2 = data_importer->image4d->getYDim();
                size3= dim3 = data_importer->image4d->getZDim();
                size4= dim4 = data_importer->image4d->getCDim();
                size5= dim5 = 1;
                if (data_importer->image4d->getTDim()>1 && data_importer->image4d->getTimePackType()==TIME_PACK_C)
                {
                    MESSAGE_ASSERT(data_importer->image4d->getCDim() >= data_importer->image4d->getTDim());

                    size4=dim4 = data_importer->image4d->getCDim()/data_importer->image4d->getTDim();
                    size5=dim5 = data_importer->image4d->getTDim();
                }
                start1 = 0;
                start2 = 0;
                start3 = 0;
                start4 = 0;
                start5 = 0;
            }
            else // image4d==0  default coordinate frame for surface
            {
                size1=dim1 = 0; //DEFAULT_DIM1;
                size2=dim2 = 0; //DEFAULT_DIM2;
                size3=dim3 = 0; //DEFAULT_DIM3;
                size4=dim4 = 0; // this make no rgbaBuf allocated
                size5=dim5 = 0; // packed time
                start1 = 0;
                start2 = 0;
                start3 = 0;
                start4 = 0;
                start5 = 0;
                //data4dp = 0; // this make no rgbaBuf allocated
            }


            //         if (b_limitedsize)
            //         {
            //              getLimitedSampleScaleBufSize(size1, size2, size3, size4, size5, sampleScale, bufSize);
            //         }
            //         else
            //         {
            bufSize[0] = size1;
            bufSize[1] = size2;
            bufSize[2] = size3;
            bufSize[3] = size4;
            bufSize[4] = size5;
            //              sampleScale[0]=sampleScale[1]=sampleScale[2]=sampleScale[3]=sampleScale[4] = 1;
            //         }

            total_rgbaBuf = rgbaBuf = 0; //(RGBA*)-1; //test whether the new sets pointer to 0 when failed
            if (size4>0)
            {
                // only RGB, first 3 channels of original image
                //total_rgbaBuf = rgbaBuf = new RGBA8[ bufSize[0] * bufSize[1] * bufSize[2] * 1 * bufSize[4] ];
                total_rgbaBuf = new RGBA8[ bufSize[0] * bufSize[1] * bufSize[2] * 1 * bufSize[4] ];
            }

            //         qDebug("   data4dp = %0p \t(start %dx%dx%d_%d_%d, size %dx%dx%d_%d_%d)", data4dp,
            //          start1,start2,start3,start4,start5,  size1,size2,size3,size4, size5);

            qDebug("   rgbaBuf = %0p \t(%dx%dx%d_%d_%d)", rgbaBuf, bufSize[0],bufSize[1],bufSize[2],bufSize[3],bufSize[4]);


            //dataViewProcBox = dataBox = BoundingBox(start1, start2, start3, start1+(size1-1), start2+(size2-1), start3+(size3-1));

            //qDebug("   data box in original image space @\t(%g %g %g)--(%g %g %g)", dataBox.x0,dataBox.y0,dataBox.z0, dataBox.x1,dataBox.y1,dataBox.z1);

        } CATCH_handler( "RayCastCanvas:setVolume" );
        if (data_importer->image4d)
        {
            Image4DProxy<Image4DSimple> img4dp( data_importer->image4d );
            img4dp.set_minmax(data_importer->p_vmin, data_importer->p_vmax);

            timer.start();
            data_importer->data4dp_to_rgba3d(img4dp,  dim5,
                                             start1, start2, start3, start4,
                                             size1, size2, size3, size4,
                                             total_rgbaBuf, bufSize);
        }
        if (dim4==1)   data_importer->rgba3d_r2gray(total_rgbaBuf, bufSize); //081103
        qDebug() << "The slow operation took" << ((float)timer.elapsed())/1000.0 << "seconds";
    }
    double p_min = data_importer->p_vmin[0];
    double p_max = data_importer->p_vmax[0];
//    long sx = data_importer->image4d->getXDim();
//    long sy = data_importer->image4d->getYDim();
//    long sz = data_importer->image4d->getZDim();
//    long sc = data_importer->image4d->getCDim(); // for gray image stacked in channel, sc is the time indeed

    if (frame4display>=size5){
        throw std::runtime_error("data to show is not gray.");
    }
    long offsets = frame4display*size1*size2*size3;
    //rgbaBuf = total_rgbaBuf + offsets;

    if (!m_raycasting_volume)
        m_raycasting_volume = new RayCastVolume(this);
    m_raycasting_volume->transfer_volume(total_rgbaBuf + offsets, p_min, p_max, size1, size2, size3, 4/*rgbaBuf contains 4 channels*/);
    update();
}
/*!
 * \brief Convert a mouse position into normalised canvas coordinates.
 * Normalized coordinates: Center of the canvas is the origin piont.
 * left-up corner is (-1,1) and right-bottom is (1,-1)
 * \param p Mouse position.
 * \return Normalised coordinates for the mouse position.
 */
QPointF RayCastCanvas::pixel_pos_to_view_pos(const QPointF& p)
{
    return QPointF(2.0 * float(p.x()) / width() - 1.0,
                   1.0 - 2.0 * float(p.y()) / height());
}

void RayCastCanvas::setLightPositionZero(){
    m_trackBall.reset2origin();
    centerShift->setX(0);
    centerShift->setY(0);
    //initializeGL();
    update();
}
/*!
 * \brief Callback for mouse movement.
 */
void RayCastCanvas::mouseMoveEvent(QMouseEvent *event)
{
    //float  test = 0;
    if (event->buttons() & Qt::LeftButton) {
        m_trackBall.move(pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        //} else {
    //    m_trackBall.release(pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
    }

}


/*!
 * \brief Callback for mouse press.
 */
void RayCastCanvas::mousePressEvent(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton) {
        m_trackBall.push(pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        m_trackBall.start();
    }
}


/*!
 * \brief Callback for mouse release.
 */
void RayCastCanvas::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_trackBall.release(pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        m_trackBall.stop();
    }
}


/*!
 * \brief Callback for mouse wheel.
 */
void RayCastCanvas::wheelEvent(QWheelEvent * event)
{
    m_distExp += event->angleDelta().y(); // for modern mouse, there may be two wheels; common one has only vertical one using angleDelta().y()
    if (m_distExp < -1800)
        m_distExp = -1800;
    if (m_distExp > 600)
        m_distExp = 600;
    update();
}
/*!
 * \brief slot: display frame t.
 */
void RayCastCanvas::setVolumeTimePoint(int t)
{
    //qDebug("V3dR_GLWidget::setVolumeTimePoint = %d", t);
    if (t<0) t = 0;
    if (t>=data_importer->image4d->getTDim()){
        t = data_importer->image4d->getTDim()-1;
    }
    this->setVolume(t);
    emit changeVolumeTimePoint(t); //need?
}
/*!
 * \brief Add a shader.
 * \param name Name for the shader.
 * \param vertex Vertex shader source file.
 * \param fragment Fragment shader source file.
 */
void RayCastCanvas::add_shader(const QString& name, const QString& vertex, const QString& fragment)
{
    m_shaders[name] = new QOpenGLShaderProgram(this);
    m_shaders[name]->addShaderFromSourceFile(QOpenGLShader::Vertex, vertex);
    m_shaders[name]->addShaderFromSourceFile(QOpenGLShader::Fragment, fragment);
    m_shaders[name]->link();
}
/*!
 * \brief slot: change the contrast.
 */
void RayCastCanvas::setContrast(int relative_contrast/*[-100:100]*/)
{
    if (relative_contrast == 0) m_gamma = 1.0;
    else if (relative_contrast > 0) m_gamma = 1.0 + relative_contrast/25.0;
    else m_gamma = 1.0 + relative_contrast * 0.008;
    //m_gamma = (relative_contrast+100.0)/40.0;
    update();
    //RayCastCanvas::raycasting(const QString& shader);
}
void RayCastCanvas::handleKeyPressEvent(QKeyEvent * e)  //090428 RZC: make public function to finally overcome the crash problem of hook MainWindow
{
    //qDebug("size: %d, %f, %d, %f", scaled_width(), ((float)scaled_width()), scaled_height(), ((float)scaled_height()));
    float stepSizeX = 20.0/((float)scaled_width());
    float stepSizeY = 20.0/((float)scaled_height());
    bool update_flag = true;
    switch (e->key())
    {
        case Qt::Key_Left: //100802: arrows key must use WITH_?_MODIFIER
            {
                //centerShift->x()-1;
                centerShift->setX(centerShift->x()-stepSizeX);
            }
            break;
        case Qt::Key_Right:
            {
                centerShift->setX(centerShift->x()+stepSizeX);
            }
            break;
        case Qt::Key_Up:
            {
                centerShift->setY(centerShift->y()+stepSizeY);
            }
            break;
        case Qt::Key_Down:
            {
                centerShift->setY(centerShift->y()-stepSizeY);
            }
            break;

            //////////////////////////////////////////////////////////////////////////////
        default:
            update_flag = false;
            //QOpenGLWidget::keyPressEvent(e);
            break;
    }
    if (update_flag) update();
}
void RayCastCanvas::handleKeyReleaseEvent(QKeyEvent * e)  //090428 RZC: make public function to finally overcome the crash problem of hook MainWindow
{
    QOpenGLWidget::keyReleaseEvent(e);
    update(); //091030: must be here for correct MarkerPos's view matrix
}
void RayCastCanvas::setBnfAxesOnOff()
{
    if(m_raycasting_volume){
        bShowBoundingBox = !bShowBoundingBox;
    }
}
/**rendering text*/
void RayCastCanvas::renderText(double x, double y, double z, QString text)
{
    int width = this->width();
    int height = this->height();

    GLdouble model[4][4], proj[4][4];
    GLint view[4];
    glGetDoublev(GL_MODELVIEW_MATRIX, &model[0][0]);
    glGetDoublev(GL_PROJECTION_MATRIX, &proj[0][0]);
    glGetIntegerv(GL_VIEWPORT, &view[0]);
    GLdouble textPosX = 0, textPosY = 0, textPosZ = 0;

    project(x, y, z,
                &model[0][0], &proj[0][0], &view[0],
                &textPosX, &textPosY, &textPosZ);

    textPosY = height - textPosY; // y is inverted

    QPainter painter(this);
    painter.setPen(Qt::yellow);
    painter.setFont(QFont("Helvetica", 8));
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    painter.drawText(textPosX, textPosY, text); // z = pointT4.z + distOverOp / 4
    painter.end();
}

inline GLint RayCastCanvas::project(GLdouble objx, GLdouble objy, GLdouble objz,
    const GLdouble model[16], const GLdouble proj[16],
    const GLint viewport[4],
    GLdouble * winx, GLdouble * winy, GLdouble * winz)
{
    GLdouble in[4], out[4];

    in[0] = objx;
    in[1] = objy;
    in[2] = objz;
    in[3] = 1.0;
    transformPoint(out, model, in);
    transformPoint(in, proj, out);

    if (in[3] == 0.0)
        return GL_FALSE;

    in[0] /= in[3];
    in[1] /= in[3];
    in[2] /= in[3];

    *winx = viewport[0] + (1 + in[0]) * viewport[2] / 2;
    *winy = viewport[1] + (1 + in[1]) * viewport[3] / 2;

    *winz = (1 + in[2]) / 2;
    return GL_TRUE;
}

inline void RayCastCanvas::transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4])
{
#define M(row,col)  m[col*4+row]
    out[0] =
        M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
    out[1] =
        M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
    out[2] =
        M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
    out[3] =
        M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
#undef M
}


void RayCastCanvas::cleanData()
{
    qDebug("   Renderer_gl1::cleanData");

    for (int i=0; i<5; i++)
    {
        //sampleScale[i]=1;
        bufSize[i]=0;
    }

    DELETE_AND_ZERO(total_rgbaBuf);
    rgbaBuf = 0;
    //DELETE_AND_ZERO(rgbaBuf_Yzx);
    //DELETE_AND_ZERO(rgbaBuf_Xzy);
}


/*** The following functions are for bounding box and axes draw ****/
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
void RayCastCanvas::draw_bbox() {
    qDebug("%d, %f \n", this->width(), this->m_raycasting_volume->get_size().x());
    float exceed_extent = 1.2;
    float td = 0.015;
    QVector3D newStartPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, 1);
    QVector3D newEndPt = m_modelViewProjectionMatrix * QVector3D(exceed_extent, -1, 1);
    glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
    {
        glColor3f(1, 0, 0);
        glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
        glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
    }
    glEnd();

    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
             XYZ(0.0f, 0.0f, td));
    glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
    {
        glColor3f(0, 1, 0);
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, exceed_extent, 1);
        glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
        glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
    }
    glEnd();
    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
             XYZ(0.0f, 0.0f, td));
    glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
    {
        glColor3f(0, 0, 1);
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, -exceed_extent);
        glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
        glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
    }
    glEnd();
    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
             XYZ(0.0f, 0.0f, td));


    // Cube 1x1x1, centered on origin
//    GLfloat vertices[] = {
//        -1.5, -0.5, -0.5, 1.0,
//        1.5, -0.5, -0.5, 1.0,
//        1.5,  0.5, -0.5, 1.0,
//        -1.5,  0.5, -0.5, 1.0,
//        -1.5, -0.5,  0.5, 1.0,
//        1.5, -0.5,  0.5, 1.0,
//        1.5,  0.5,  0.5, 1.0,
//        -1.5,  0.5,  0.5, 1.0,
//    };
//    GLuint vbo_vertices;
//    glGenBuffers(1, &vbo_vertices);
//    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
//    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
//    glBindBuffer(GL_ARRAY_BUFFER, 0);

//    GLushort elements[] = {
//        0, 1, 2, 3,
//        4, 5, 6, 7,
//        0, 4, 1, 5, 2, 6, 3, 7
//    };
//    GLuint ibo_elements;
//    glGenBuffers(1, &ibo_elements);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);


//    //  GLfloat
//    //    min_x, max_x,
//    //    min_y, max_y,
//    //    min_z, max_z;
//    //  min_x = max_x = mesh->vertices[0].x;
//    //  min_y = max_y = mesh->vertices[0].y;
//    //  min_z = max_z = mesh->vertices[0].z;
//    //  for (int i = 0; i < mesh->vertices.size(); i++) {
//    //    if (mesh->vertices[i].x < min_x) min_x = mesh->vertices[i].x;
//    //    if (mesh->vertices[i].x > max_x) max_x = mesh->vertices[i].x;
//    //    if (mesh->vertices[i].y < min_y) min_y = mesh->vertices[i].y;
//    //    if (mesh->vertices[i].y > max_y) max_y = mesh->vertices[i].y;
//    //    if (mesh->vertices[i].z < min_z) min_z = mesh->vertices[i].z;
//    //    if (mesh->vertices[i].z > max_z) max_z = mesh->vertices[i].z;
//    //  }
//    //  glm::vec3 size = glm::vec3(max_x-min_x, max_y-min_y, max_z-min_z);
//    //  glm::vec3 center = glm::vec3((min_x+max_x)/2, (min_y+max_y)/2, (min_z+max_z)/2);
//    //  glm::mat4 transform = glm::translate(glm::mat4(1), center) * glm::scale(glm::mat4(1), size);

//    //  /* Apply object's transformation matrix */
//    //  glm::mat4 m = mesh->object2world * transform;
//    //  glUniformMatrix4fv(uniform_m, 1, GL_FALSE, glm::value_ptr(m));
//    GLint attribute_v_coord = -1;
//    glBindBuffer(GL_ARRAY_BUFFER, vbo_vertices);
//    glEnableVertexAttribArray(attribute_v_coord);
//    glVertexAttribPointer(
//                attribute_v_coord,  // attribute
//                4,                  // number of elements per vertex, here (x,y,z,w)
//                GL_FLOAT,           // the type of each element
//                GL_FALSE,           // take our values as-is
//                0,                  // no extra data between each position
//                0                   // offset of first element
//                );

//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_elements);
//    glLineWidth(3); glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, 0);
//    glLineWidth(3); glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_SHORT, (GLvoid*)(4*sizeof(GLushort)));
//    glLineWidth(3); glDrawElements(GL_LINES, 8, GL_UNSIGNED_SHORT, (GLvoid*)(8*sizeof(GLushort)));
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

//    glDisableVertexAttribArray(attribute_v_coord);
//    glBindBuffer(GL_ARRAY_BUFFER, 0);

//    glDeleteBuffers(1, &vbo_vertices);
//    glDeleteBuffers(1, &ibo_elements);
}


void RayCastCanvas::setBoundingBoxSpace(BoundingBox BB)
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



void RayCastCanvas::drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth, float AlineWidth)
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
        glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
        {
//            glColor3f(1, 0, 0);		box_quads( BoundingBox(A0, XYZ(A1.x, A0.y+ld, A0.z+ld)) );
//            glColor3f(0, 1, 0);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A1.y, A0.z+ld)) );
//            glColor3f(0, 0, 1);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A0.y+ld, A1.z)) );
            //glColor3f(1, 0, 0);
            glVertex3f(0, 0, 0); glVertex3f(10, 0, 0);
            //glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 10, 0);
            //glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 10);
        }
        glEnd();

//        glColor3f(1, 0, 0);		drawString(A1.x+td, A0.y, A0.z, "X", 1, 0);
//        glColor3f(0, 1, 0);		drawString(A0.x, A1.y+td, A0.z, "Y", 1, 0);
//        glColor3f(0, 0, 1);		drawString(A0.x, A0.y, A1.z+td, "Z", 1, 0);
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
void RayCastCanvas::drawString(float x, float y, float z, const char* text, int shadow, int fontsize)
{
    //if (! this)  return;

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
            QPainter painter(this);
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
