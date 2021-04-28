#include <cstdlib>
#include <ctime>
#include <iostream>
#include <vector>
#include <QCoreApplication>
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
    m_background(QColor(0, 0, 0)),
    m_raycasting_volume {nullptr},
    m_active_mode("MIP")
{
    // Register the rendering modes here, so they are available to the UI when it is initialised
    m_modes["Isosurface"] = [&]() { RayCastCanvas::raycasting("Isosurface"); };
    m_modes["Alpha blending"] = [&]() { RayCastCanvas::raycasting("Alpha blending"); };
    m_modes["Alpha blending rgba"] = [&]() { RayCastCanvas::raycasting("Alpha blending rgba"); };
    m_modes["MIP"] = [&]() { RayCastCanvas::raycasting("MIP"); };

    m_init_mode = m_active_mode;
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
    //glEnable(GL_MULTISAMPLE);
    if(hasMouseTracking())
        setMouseTracking(false);

    m_raycasting_volume = new RayCastVolume(this);

    add_shader("Isosurface", ":/shaders/isosurface.vert", ":/shaders/isosurface.frag");
    add_shader("Alpha blending", ":/shaders/alpha_blending.vert", ":/shaders/alpha_blending.frag");
    add_shader("Alpha blending rgba", ":/shaders/alpha_blending.vert", ":/shaders/alpha_blending_rgba.frag");
    add_shader("MIP", ":/shaders/maximum_intensity_projection.vert", ":/shaders/maximum_intensity_projection.frag");

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
    m_raycasting_volume->create_noise();
}

void RayCastCanvas::drawInstructions(QPainter *painter)
{
    QString text = tr("Click and drag with the left mouse button "
                      "to rotate the data.");
    QFontMetrics metrics = QFontMetrics(font());
    int border = qMax(4, metrics.leading());

    QRect rect = metrics.boundingRect(0, 0, width() - 2*border, int(height()*0.125),
                                      Qt::AlignCenter | Qt::TextWordWrap, text);
    painter->setRenderHint(QPainter::TextAntialiasing);
    painter->fillRect(QRect(0, 0, width(), rect.height() + 2*border),
                     QColor(0, 0, 0, 127));
    painter->setPen(Qt::white);
    painter->fillRect(QRect(0, 0, width(), rect.height() + 2*border),
                      QColor(0, 0, 0, 127));
    painter->drawText((width() - rect.width())/2, border,
                      rect.width(), rect.height(),
                      Qt::AlignCenter | Qt::TextWordWrap, text);

}
/**
 * @brief RayCastCanvas::drawLine at a specific location related to the rendered volume
 * @param painter
 * @param c
 * @param p0: start location in camera view, which x,y,z in [-1, 1]
 * @param p1: end location in camera view, which x,y,z in [-1, 1]
 * @param lineWidth
 */
void RayCastCanvas::drawLine(QPainter *painter, QColor c, QPointF p0, QPointF p1 , int lineWidth)
{
    painter->setPen(QPen(c, lineWidth));
    QPointF p_start = view_pos_to_pixel_pos(p0);
    QPointF p_end = view_pos_to_pixel_pos(p1);
    painter->drawLine(p_start, p_end);

}
/**
 * @brief RayCastCanvas::drawText at a specific location related to the rendered volume
 * @param painter
 * @param c
 * @param p: location in camera view, which x,y,z in [-1, 1]
 * @param text
 */
void RayCastCanvas::drawText(QPainter *painter, QColor c, QPointF p, QString text){
    //QVector3D newStartPt = m_modelViewProjectionMatrix * QVector3D(1.3, 1, 1);
    painter->setPen(c);
    painter->setFont(QFont("Arial", 10)); //font not sure large enough or not
    QPointF pAtScreen = view_pos_to_pixel_pos(p);
    painter->drawText(pAtScreen.x(), pAtScreen.y(), text);

}
/*!
 * \brief Paint a frame on the canvas.
 */
void RayCastCanvas::paintGL()
{
    if (!m_raycasting_volume && data_importer) {
        //makeCurrent();
        m_raycasting_volume = new RayCastVolume(this);
    }

    if (m_raycasting_volume){
        m_raycasting_volume->create_noise();
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
        //glPushMatrix();
        m_modes[m_active_mode]();
        //glPopMatrix();
    }
    if (bShowAxes)
    {
        glPushMatrix();
        draw_axes();
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
    if (m_raycasting_volume){
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
            m_shaders[shader]->setUniformValue("consider_transparency", m_consider_transparency /*0*/);
            m_shaders[shader]->setUniformValue("min_valid_intensity", m_min_valid_intensity /*0*/);
            m_shaders[shader]->setUniformValue("leftup_xyz", m_leftup_xyz /*0*/);
            m_shaders[shader]->setUniformValue("rightbottom_xyz", m_rightbottom_xyz /*0*/);

            glClearColor(m_background.redF(), m_background.greenF(), m_background.blueF(), m_background.alphaF());
            glClear(GL_COLOR_BUFFER_BIT);

            m_raycasting_volume->paint();
        }
        m_shaders[shader]->release();
    }
}
/*!
 * \brief import the data and set it as color RGBA data if needed
 */
void RayCastCanvas::setVolume(long frame4display) {

    if(show_track_result){ // display RGBA tracking results
        if(rgb_frame.empty()){
            throw std::runtime_error("rgb results has not been initialized.");
        }
        this->setMode("Alpha blending rgba");//
        //m_gamma = 1;
        m_raycasting_volume->transfer_volume((unsigned char *)rgb_frame.data, 0, 255, rgb_frame.size[1],
                rgb_frame.size[0], rgb_frame.size[2], 4);
    }else{
        if (!data_importer)
        {
            throw std::runtime_error("data_importer has not been initialized.");
        }
        if (!data_importer->p_vmin){// if max min value not defined
            data_importer->updateminmaxvalues();
        }
        curr_timePoint_in_canvas = frame4display; // default is 0
        try
        {
            if (data_importer->image4d && data_importer->image4d->getCDim()>0)
            {
                bufSize[0] = data_importer->image4d->getXDim();
                bufSize[1] = data_importer->image4d->getYDim();
                bufSize[2] = data_importer->image4d->getZDim();
                bufSize[3] = data_importer->image4d->getCDim();
                bufSize[4] = 1;
            }
            if (data_importer->image4d->getTDim()>1 && data_importer->image4d->getTimePackType()==TIME_PACK_C)
            {
                MESSAGE_ASSERT(data_importer->image4d->getCDim() >= data_importer->image4d->getTDim());

                bufSize[3] = data_importer->image4d->getCDim()/data_importer->image4d->getTDim();
                bufSize[4] = data_importer->image4d->getTDim();
            }
        } CATCH_handler( "RayCastCanvas:setVolume" );
        if (!total_rgbaBuf && (bufSize[3] == 3 || flag_rgba_display)) // if buffer has not been initialized
        {
            QElapsedTimer timer;
            // now the correct data size xyzct are saved in bufSize
            total_rgbaBuf = rgbaBuf = 0; //(RGBA*)-1; //test whether the new sets pointer to 0 when failed
            if (bufSize[3]>0)
            {
                // only RGB, first 3 channels of original image
                total_rgbaBuf = new RGBA8[ bufSize[0] * bufSize[1] * bufSize[2] * 1 * bufSize[4]];
            }

            if (data_importer->image4d)
            {
                Image4DProxy<Image4DSimple> img4dp( data_importer->image4d );
                img4dp.set_minmax(data_importer->p_vmin, data_importer->p_vmax);

                timer.start();
                data_importer->data4dp_to_rgba3d(img4dp,  bufSize[4],
                        0, 0, 0, 0, // the starting location to make the transfer
                        bufSize[0], bufSize[1], bufSize[2], bufSize[3],
                        total_rgbaBuf, bufSize);
            }
            if (data_importer->image4d->getCDim()==1)   data_importer->rgba3d_r2gray(total_rgbaBuf, bufSize); //081103
            qDebug() << "The slow operation took" << ((float)timer.elapsed())/1000.0 << "seconds";
        }

        double p_min = data_importer->p_vmin[0];
        double p_max = data_importer->p_vmax[0];

        if (frame4display>=bufSize[4]){
            throw std::runtime_error("data to show is not gray.");
        }
        long offsets = frame4display*bufSize[0]*bufSize[1]*bufSize[2];
        if (total_rgbaBuf){ // display the data using 4 channels

            m_raycasting_volume->transfer_volume(total_rgbaBuf + offsets,
                                                 p_min, p_max, bufSize[0],
                    bufSize[1], bufSize[2], 4/*rgbaBuf contains 4 channels*/);
            this->setMode("Alpha blending rgba");//
//            if(frame4display == 1){
//                resetMode(); // reset the renderering to gray-scale
//            }else{
//                this->setMode("Alpha blending rgba");//
//            }
        }else{
            offsets *= data_importer->image4d->getUnitBytes(); // if 16 bit, one pixel occupies two chars.
            m_raycasting_volume->transfer_volume(data_importer->image4d->getRawData() + offsets,
                                                 p_min, p_max, bufSize[0],
                    bufSize[1], bufSize[2], 1);
            resetMode(); // reset the renderering to gray-scale
        }

    }
    if(!m_gamma_init){
        //double max_exist_intensity = data_importer->p_vmax[0];
        double max_intensity;
        if (data_importer->image4d->getDatatype() == V3D_UINT8){
            max_intensity = 255;
        }else if(data_importer->image4d->getDatatype() == V3D_UINT16){
            max_intensity = 65535;
        }else{
            max_intensity = 0;
            qDebug("unsupported image type");
        }
        if(data_importer->p_vmax[0] > 10 && data_importer->p_vmax[0] < max_intensity/4){
            m_gamma = log(data_importer->p_vmax[0]/max_intensity) / log(0.25);
        }else{
            m_gamma = 1.0;
        }
        m_gamma_init = true;
        m_gamma0 = m_gamma;
    }
    update();


}


/**
 * @brief RayCastCanvas::setVolumeWithMask
 * @param frame4display
 * @param mask
 * @param datatype
 */
void RayCastCanvas::setVolumeWithMask(long frame4display, unsigned char* mask) {
    //datatype == 1: unsigned char
    //datatype == 2: int 32
    if (!data_importer)
    {
        throw std::runtime_error("data_importer has not been initialized.");
    }
    if (!data_importer->p_vmin){// if max min value not defined
        data_importer->updateminmaxvalues();
    }
    size_t numVox = bufSize[0] * bufSize[1] * bufSize[2];
    RGBA8 *imWithMask = new RGBA8[ numVox ];
    double p_min = data_importer->p_vmin[frame4display];
    double p_max = data_importer->p_vmax[frame4display];

    // scale the image to 0, 1
    long offsets = frame4display*bufSize[0]*bufSize[1]*bufSize[2];
    offsets *= data_importer->image4d->getUnitBytes();
    void *start = data_importer->image4d->getRawData() + offsets;
    float vox_val = 0;
    for (size_t i = 0; i < numVox; i++){
        if (mask[i] > 0){
            imWithMask[i].r = 150;
        }else{
            if (data_importer->image4d->getUnitBytes() == 1){
                vox_val = (float)*((v3d_uint8 *)start + i);
            }else{
                vox_val = (float)*((v3d_uint16 *)start + i);
            }
        }
        vox_val = 2*255.0*(vox_val-p_min)/(p_max-p_min); // times 2 to make it clearer
        imWithMask[i].g = (v3d_uint8)vox_val;
        imWithMask[i].b = 0;
        imWithMask[i].a = (v3d_uint8)(0.0f + imWithMask[i].r + imWithMask[i].g + + imWithMask[i].b);
    }

    // display the data using 4 channels
    m_raycasting_volume->transfer_volume(imWithMask,
                                         p_min, p_max, bufSize[0],
            bufSize[1], bufSize[2], 4/*rgba contains 4 channels*/);
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
QPointF RayCastCanvas::view_pos_to_pixel_pos(const QPointF& p)
{
    return QPointF(float(p.x() + 1) * width() / 2.0,
                   height() * float(1.0 - p.y()) / 2.0);
}
/**
 * @brief RayCastCanvas::setLightPositionZero
 * Reset the camera position and light position
 */
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
void RayCastCanvas::setVolumeTimePoint(int t)//,
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
    //if (relative_contrast == 0) m_gamma = m_gamma;
    //else
    if (relative_contrast > 0) m_gamma = m_gamma0 + relative_contrast/25.0;
    else m_gamma = m_gamma0 + relative_contrast * 0.008;
    //m_gamma = (relative_contrast+100.0)/40.0;
    update();
    //RayCastCanvas::raycasting(const QString& shader);
}
void RayCastCanvas::setThreshold(int intensity_threshold){
//    double p_min = data_importer->p_vmin[0];
//    double p_max = data_importer->p_vmax[0];
    double max_intensity;
    if (data_importer->image4d->getDatatype() == V3D_UINT8){
        max_intensity = 255;
    }else if(data_importer->image4d->getDatatype() == V3D_UINT16){
        max_intensity = 65535;
    }else{
        max_intensity = 0;
        qDebug("unsupported image type");
    }
    if(intensity_threshold == 0){
        m_min_valid_intensity = 0;
    }else if(intensity_threshold == 100){
        m_min_valid_intensity = 1;
    }else{
        m_min_valid_intensity = pow(intensity_threshold,3) / 1000000.0;
        m_min_valid_intensity *= data_importer->p_vmax[bufSize[4]-1] / max_intensity ;
    }
    update();

}

void RayCastCanvas::setRangeXMIN(int xmin){
    m_leftup_xyz.setX(xmin / 100.0);
    update();
}
void RayCastCanvas::setRangeXMAX(int xmax){
    m_rightbottom_xyz.setX(xmax / 100.0);
    update();
}

void RayCastCanvas::setRangeYMIN(int ymin){
    m_leftup_xyz.setY(ymin / 100.0);
    update();
}

void RayCastCanvas::setRangeYMAX(int ymax){
    m_rightbottom_xyz.setY(ymax / 100.0);
    update();
}

void RayCastCanvas::setRangeZMIN(int zmin){ //z-direction is different
    m_rightbottom_xyz.setZ(1 - zmin / 100.0);
    update();
}

void RayCastCanvas::setRangeZMAX(int zmax){//z-direction is different
    m_leftup_xyz.setZ(1 - zmax / 100.0);
    update();
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
        bShowAxes = !bShowAxes;
    }
    //emit changeBnfAxesOnOff(bShowAxes);
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
/**
 * @brief draw_tri: draw a filled triangle (src from vaa3d), (deprecated)
 * @param P1
 * @param P2
 * @param P3
 * @param offst
 */
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
/**
 * @brief draw_axes, add axes in the figure with text X, Y, and Z
 */
void RayCastCanvas::draw_axes() {
        float exceed_extent = 1.2;
        //qDebug ("drawInstructions has been called \n"); //this->paintText();
        QPainter painter(this);
        QVector3D originPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, 1);
        // x-axis
        QVector3D newEndPt = m_modelViewProjectionMatrix * QVector3D(exceed_extent, -1, 1);
        QColor c = QColor(255,0,0);
        this->drawLine(&painter, c, QPointF(originPt.x(), originPt.y()),
                       QPointF(newEndPt.x(), newEndPt.y()), 2);
        this->drawText(&painter, c, QPointF(newEndPt.x(), newEndPt.y()), "X");
        // y-axis
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, exceed_extent, 1);
        c = QColor(0,255,0);
        this->drawLine(&painter, c, QPointF(originPt.x(), originPt.y()),
                       QPointF(newEndPt.x(), newEndPt.y()), 2);
        this->drawText(&painter, c, QPointF(newEndPt.x(), newEndPt.y()), "Y");
        // z-axis
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, -exceed_extent);
        c = QColor(0,0,255);
        this->drawLine(&painter, c, QPointF(originPt.x(), originPt.y()),
                       QPointF(newEndPt.x(), newEndPt.y()), 2);
        this->drawText(&painter, c, QPointF(newEndPt.x(), newEndPt.y()), "Z");
        // draw instructions
        //this->drawInstructions(&painter);
        painter.end();

        if (false){// deprecated way to draw lines
            qDebug("%d, %f \n", this->width(), this->m_raycasting_volume->get_size().x());
            float exceed_extent = 1.2;
            //float td = 0.015;
            QVector3D newStartPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, 1);
            QVector3D newEndPt = m_modelViewProjectionMatrix * QVector3D(exceed_extent, -1, 1);
            glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
            {
                glColor3f(1, 0, 0);//glLineWidth(3);
                glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
                glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
            }
            glEnd();
            //    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
            //             XYZ(0.0f, 0.0f, td));
            glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
            {
                glColor3f(0, 1, 0);//glLineWidth(3);
                newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, exceed_extent, 1);
                glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
                glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
            }
            glEnd();
            //    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
            //             XYZ(0.0f, 0.0f, td));
            glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
            {
                glColor3f(0, 0, 1); //glLineWidth(3);
                newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, -1, -exceed_extent);
                glVertex3f(newStartPt.x(), newStartPt.y(), newStartPt.z());
                glVertex3f(newEndPt.x(), newEndPt.y(), newEndPt.z());
            }
            glEnd();
            //    draw_tri(XYZ(newEndPt.x()-td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+3.0f*td, newEndPt.y()+td, newEndPt.z()),
            //             XYZ(newEndPt.x()+td, newEndPt.y()+3.0f*td, newEndPt.z()),
            //             XYZ(0.0f, 0.0f, td));
        }
}

/**
 * @brief RayCastCanvas::renderText-> the following 3 functions
 *  are a way to render text using QPainter. keep it as a back up method.
 * @param x
 * @param y
 * @param z
 * @param text
 */
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

/** this is the Vaa3d's idea of adding axes. It fails on our framework. Use Qpainter instead.*/
//void RayCastCanvas::setBoundingBoxSpace(BoundingBox BB)
//{
//    float DX = BB.Dx();
//    float DY = BB.Dy();
//    float DZ = BB.Dz();
//    float maxD = BB.Dmax();

//    double s[3];
//    s[0] = 1/maxD *2;
//    s[1] = 1/maxD *2;
//    s[2] = 1/maxD *2;
//    double t[3];
//    t[0] = -BB.x0 -DX /2;
//    t[1] = -BB.y0 -DY /2;
//    t[2] = -BB.z0 -DZ /2;

//    // from boundingBox space ==> fit in [-1, +1]^3
//    glScaled(s[0], s[1], s[2]);
//    glTranslated(t[0], t[1], t[2]);
//}

//inline void box_quads(const BoundingBox & BB)
//{
//#define BB_VERTEX(xi,yi,zi)  glVertex3d(BB.x##xi, BB.y##yi, BB.z##zi)

//    BB_VERTEX(0, 0, 0);	BB_VERTEX(0, 1, 0);	BB_VERTEX(1, 1, 0); BB_VERTEX(1, 0, 0); //z=0
//    BB_VERTEX(0, 0, 1);	BB_VERTEX(0, 1, 1);	BB_VERTEX(1, 1, 1); BB_VERTEX(1, 0, 1); //z=1

//    BB_VERTEX(0, 0, 0);	BB_VERTEX(1, 0, 0);	BB_VERTEX(1, 0, 1); BB_VERTEX(0, 0, 1); //y=0
//    BB_VERTEX(0, 1, 0);	BB_VERTEX(1, 1, 0);	BB_VERTEX(1, 1, 1); BB_VERTEX(0, 1, 1); //y=1

//    BB_VERTEX(0, 0, 0);	BB_VERTEX(0, 0, 1);	BB_VERTEX(0, 1, 1); BB_VERTEX(0, 1, 0); //x=0
//    BB_VERTEX(1, 0, 0);	BB_VERTEX(1, 0, 1);	BB_VERTEX(1, 1, 1); BB_VERTEX(1, 1, 0); //x=1

//}



//void RayCastCanvas::drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth, float AlineWidth)
//{
//    glPushAttrib(GL_LINE_BIT | GL_POLYGON_BIT);
//            //| GL_DEPTH_BUFFER_BIT);
//    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//    glEnable(GL_POLYGON_OFFSET_LINE);

////	glPolygonOffset(0, -1); // deal z-fighting, 081120
////	glDepthFunc(GL_LEQUAL);
//    if (posXTranslateBB != 0) delete posXTranslateBB;
//    if (negXTranslateBB != 0) delete negXTranslateBB;
//    if (posYTranslateBB != 0) delete posYTranslateBB;
//    if (negYTranslateBB != 0) delete negYTranslateBB;
//    posXTranslateBB=0, negXTranslateBB=0,
//            posYTranslateBB=0, negYTranslateBB=0;
//    // an indicator of coordinate direction
//    if (bShowAxes && AlineWidth>0)
//    {
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

//        float D = (BB.Dmax());
//        float ld = D*0.0001; //1e-4 is best
//        float td = D*0.015;
//        XYZ A0 = BB.Vabsmin();
//        XYZ A1 = BB.V1() + D*0.05;

//        glPolygonOffset(-0.002, -2); //(-0.002, -2) for good z-fighting with bounding box, 081120,100823

//        glLineWidth(AlineWidth); // work only before glBegin(), by RZC 080827
//        //glBegin(GL_QUADS);
//        glBegin(GL_LINES); // glPolygonOffset do NOT  influence GL_LINES
//        {
////            glColor3f(1, 0, 0);		box_quads( BoundingBox(A0, XYZ(A1.x, A0.y+ld, A0.z+ld)) );
////            glColor3f(0, 1, 0);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A1.y, A0.z+ld)) );
////            glColor3f(0, 0, 1);		box_quads( BoundingBox(A0, XYZ(A0.x+ld, A0.y+ld, A1.z)) );
//            //glColor3f(1, 0, 0);
//            glVertex3f(0, 0, 0); glVertex3f(10, 0, 0);
//            //glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, 10, 0);
//            //glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, 10);
//        }
//        glEnd();

////        glColor3f(1, 0, 0);		drawString(A1.x+td, A0.y, A0.z, "X", 1, 0);
////        glColor3f(0, 1, 0);		drawString(A0.x, A1.y+td, A0.z, "Y", 1, 0);
////        glColor3f(0, 0, 1);		drawString(A0.x, A0.y, A1.z+td, "Z", 1, 0);
//    }

//    if (bShowBoundingBox && BlineWidth>0)
//    {
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

//        glPolygonOffset(0, -1); // deal z-fighting with volume, 081120

//        glLineWidth(BlineWidth); // work only before glBegin(), by RZC 080827
//        glBegin(GL_QUADS);
//        //glBegin(GL_LINES);
//        {
//            glColor3fv(color_line.c);	box_quads(BB);
//        }
//        glEnd();
//    }

//    glPopAttrib();
//}
//void RayCastCanvas::drawString(float x, float y, float z, const char* text, int shadow, int fontsize)
//{
//    //if (! this)  return;

//    if (shadow)
//    {
//        glPushAttrib(GL_DEPTH_BUFFER_BIT);
//        glPushAttrib(GL_CURRENT_BIT);
//            //glColor3ub(50,50,50);
//            //glColor3ub(200,200,200);

//            // CMB MSVC debugger with Qt 4.7 triggers assert if font weight > 99
//            // QFont f;  f.setPointSize(f.pointSize()+1); f.setWeight(f.weight()+200);
//            QFont f;  f.setPointSize(f.pointSize()+1); f.setWeight(99);
////#if defined(USE_Qt5)
////#else
////            ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text), f);
////#endif
//            QPainter painter(this);
//            painter.setPen(QColor(50,50,50));
//            painter.setFont(QFont("Arial", 16));
//            painter.drawText(QPoint(int(200),int(200)), QString(text));
//            painter.end();
//        glPopAttrib();
//        glDepthFunc(GL_LEQUAL);
//    }

//    QFont f1;  f1.setPointSize((fontsize>0)?fontsize:30); f1.setWeight(99);
////    if (fontsize>0)
////#if defined(USE_Qt5)
////#else
////        ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text), f1);
////    else
////        ((QOpenGLWidget_proxy*)widget)->renderText(x,y,z, QString(text));
////#endif


//    if (shadow)
//    {
//        glPopAttrib();
//    }
//}
