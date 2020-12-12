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

        glClearColor(m_background.redF(), m_background.greenF(), m_background.blueF(), m_background.alphaF());
        glClear(GL_COLOR_BUFFER_BIT);

        m_raycasting_volume->paint();
    }
    m_shaders[shader]->release();
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


