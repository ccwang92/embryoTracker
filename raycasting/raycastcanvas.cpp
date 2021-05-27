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
    m_background(QColor(58, 90, 100)),
    m_raycasting_volume {nullptr},
    m_active_mode("MIP")
{
    // Register the rendering modes here, so they are available to the UI when it is initialised
    m_modes["Isosurface"] = [&]() { RayCastCanvas::raycasting("Isosurface"); };
    m_modes["Alpha blending"] = [&]() { RayCastCanvas::raycasting("Alpha blending"); };
    m_modes["Alpha blending rgba"] = [&]() { RayCastCanvas::raycasting("Alpha blending rgba"); };
    m_modes["MIP"] = [&]() { RayCastCanvas::raycasting("MIP"); };
    m_init_mode = m_active_mode;
    curr_timePoint_in_canvas = 0;
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

    if (total_rgbaBuf){
        delete total_rgbaBuf;
    }
    if(rgbaBuf){
        delete rgbaBuf;
    }
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
 * @brief RayCastCanvas::draw a point at a specific location related to the rendered volume
 * @param painter
 * @param c
 * @param p0: location in camera view, which x,y,z in [-1, 1]
 * @param ptWidth
 */
void RayCastCanvas::drawPoint(QPainter *painter, QColor c, QPointF p0, int ptWidth)
{
    painter->setPen(QPen(c, ptWidth));
    QPointF p_loc = view_pos_to_canvas_pixel_pos(p0);
    painter->drawPoint(p_loc);
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
    QPointF p_start = view_pos_to_canvas_pixel_pos(p0);
    QPointF p_end = view_pos_to_canvas_pixel_pos(p1);
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
    QPointF pAtScreen = view_pos_to_canvas_pixel_pos(p);
    painter->drawText(pAtScreen.x(), pAtScreen.y(), text);
}

//QVector3D RayCastCanvas::getWorldCoordinates(float mouseX, float mouseY)
//{
//    QMatrix4x4 viewMatrix ;
//    viewMatrix.lookAt(m_camera.getPosition(), m_camera.getTarget(), m_camera.getUp());

//    QMatrix4x4 projection;
//    projection.perspective(70.0, width() / height(), 0.1, 120.0);

//    QVector3D Z(0, 0, 0); // instead of 0 for x and y i need worldPosition.x() and worldPosition.y() ....
//    Z = Z.project(viewMatrix, projection, QRect(0, 0, width(), height()));

//    QVector3D worldPosition = QVector3D(x, height() - y, Z.z()).unproject(viewMatrix, projection, QRect(0, 0, width(), height()));
//    qDebug() << worldPosition;


//    QMatrix4x4 modelViewMatrix  = m_viewMatrix * m_modelMatrix;
//    QMatrix4x4 modelViewProject = m_projectionMatrix * modelViewMatrix;
//    QMatrix4x4 inverted         = viewportMatrix * modelViewProject;

//    inverted = inverted.inverted();

//    float posZ;
//    float posY = viewportSize.y() - mouseY - 1.0f;

//    world->getGLFunctions()->glReadPixels(MathHelper::toInt(mouseX), MathHelper::toInt(posY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &posZ);

//    QVector4D clickedPointOnScreen(mouseX, posY, 2.0f * posZ - 1.0f, 1.0f);
//    QVector4D clickedPointIn3DOrgn = inverted * clickedPointOnScreen;

//    clickedPointIn3DOrgn /= clickedPointIn3DOrgn.w();

//    return clickedPointIn3DOrgn.toVector3DAffine();
//}

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

        m_projectionMatrix.setToIdentity();
        m_projectionMatrix.perspective(m_fov, (float)scaled_width()/scaled_height(), 0.1f, 100.0f);
        // m_raycasting_volume->modelMatrix() saves the scale and shifting information
        m_modelMatrix = m_raycasting_volume->modelMatrix();
        m_modelViewProjectionMatrix = m_projectionMatrix * m_viewMatrix * m_modelMatrix;

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
    if(bShowMarkers){
        glPushMatrix();
        draw_makers();
        glPopMatrix();
    }

//    QPainter p(this);
//    p.setPen(Qt::red);
//    p.drawLine(rect().topLeft(), rect().bottomRight());
    if(bShowTrackResult && !bWait4RightClickOnCell){
        // we first check if rgb_frame is empty, if so, traces are already overlaid
        if ( rgb_frame.empty() && !traces.empty()){
            glPushMatrix();
            draw_traces();
            glPopMatrix();
        }
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
    if(bShowTrackResult && !rgb_frame.empty()){ // display RGBA tracking results
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
        if(data_importer){
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
 * \brief Convert a mouse position into normalised canvas coordinates (view pos).
 * Normalized coordinates: Center of the canvas is the origin piont.
 * left-up corner is (-1,1) and right-bottom is (1,-1)
 * \param p Mouse position.
 * \return Normalised coordinates for the mouse position.
 */
QPointF RayCastCanvas::canvas_pixel_pos_to_view_pos(const QPointF& p)
{
    //qDebug()<<width() << " " << scaled_width();
    //qDebug()<<height()<< " " << scaled_height();
    return QPointF(2.0 * float(p.x()) / width() - 1.0,
                   1.0 - 2.0 * float(p.y()) / height());
}
QPointF RayCastCanvas::view_pos_to_canvas_pixel_pos(const QPointF& p)
{
    return QPointF(float(p.x() + 1) * width() / 2.0,
                   height() * float(1.0 - p.y()) / 2.0);
}
/**
 * \brief Convert a voxel position in volume into normalised canvas coordinates.
 * Normalized coordinates: Center of the canvas is the origin piont.
 * left-up corner is (-1, 1, 1) and right-bottom is (1, -1, -1)
 * \param p Mouse position.
 * \return Normalised coordinates for the mouse position.
 */
QVector3D RayCastCanvas::volume_pixel_pos_to_view_pos(const QVector3D& p){
    return QVector3D(2.0 * float(p.x()) / bufSize[0] - 1.0,
                     2.0 * float(p.y()) / bufSize[1] - 1.0, //1.0 - 2.0 * float(p.y()) /  bufSize[1],
                     2.0 * float(p.z()) / bufSize[2] - 1.0);
}
QVector3D RayCastCanvas::view_pos_to_volume_pixel_pos(const QVector3D& p){
    return QVector3D(float(p.x() + 1) * bufSize[0] / 2.0,
                     float(p.y() + 1) * bufSize[1] / 2.0,//bufSize[1] * float(1.0 - p.y()) / 2.0,
                     float(p.z() + 1) * bufSize[2] / 2.0);
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
    /*! NOTE event->buttons() and event->button() are different*/
    if (event->buttons() & Qt::LeftButton) {
        m_trackBall.move(canvas_pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        //} else {
    //    m_trackBall.release(canvas_pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
    }

}


/*!
 * \brief Callback for mouse press.
 */
void RayCastCanvas::mousePressEvent(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton) {
        m_trackBall.push(canvas_pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        m_trackBall.start();
    }
}


/*!
 * \brief Callback for mouse release.
 */
void RayCastCanvas::mouseReleaseEvent(QMouseEvent *event)
{
    //qDebug() << event->button();
    if (event->button() & Qt::LeftButton) { // event->button() & Qt::LeftButton <> event->button() == Qt::LeftButton
        m_trackBall.release(canvas_pixel_pos_to_view_pos(event->pos()), m_scene_trackBall.rotation().conjugated());
        update();
        m_trackBall.stop();
    }
    if(event->button() & Qt::RightButton){
//        MarkerPos pos;
//        pos.canvas_pos = event->pos();//QPointF(tmp.x(), tmp.y());
//        pos.canvas_width = width();
//        pos.canvas_height = height();
//        pos.time_point = curr_timePoint_in_canvas;
//        pos.ModelViewProjectionMatrix = m_modelViewProjectionMatrix;
//        pos.ModelMatrix = m_modelMatrix;
//        pos.ViewMatrix = m_viewMatrix;
//        pos.ProjectionMatrix = m_projectionMatrix;
//        pos.NormalMatrix = m_normalMatrix;
//        pos.drawn = false;
//        markers.push_back(pos);
//        update();
        process_right_button_hit(event);
    }
}
/**
 * @brief userClick2CellTrace: convert user's right click to the coordinate in the volume
 * @param event
 * @return
 */
void RayCastCanvas::userClick2volumeLoc(QMouseEvent *event, QVector3D &nearEnd, QVector3D &farEnd){
    QVector3D Z(0, 0, 0); // instead of 0 for x and y i need worldPosition.x() and worldPosition.y() ....
    Z = Z.project(m_viewMatrix*m_modelMatrix,
                  m_projectionMatrix, QRect(0, 0, width(), height()));
    QVector3D loc0 = QVector3D(event->pos().x(), height() - event->pos().y(),
                      Z.z()).unproject(m_viewMatrix*m_modelMatrix, m_projectionMatrix, QRect(0, 0, width(), height()));
    qDebug() << loc0;

    Z = QVector3D(0, 0, 1);
    Z = Z.project(m_viewMatrix*m_modelMatrix,
                  m_projectionMatrix, QRect(0, 0, width(), height()));
    QVector3D loc1 =QVector3D(event->pos().x(), height() - event->pos().y(),
                              Z.z()).unproject(m_viewMatrix*m_modelMatrix, m_projectionMatrix, QRect(0, 0, width(), height()));
    qDebug() << loc1;
    // way 1: get the exact coordinate of worldPosition, e.g.: (300,200,10)
    //QVector3D curPt = worldPosition.project(m_viewMatrix*m_modelMatrix, m_projectionMatrix, QRect(0, 0, width(), height()));
    // way 2: get the normalized location of worldPosition, e.g.: (0.5, 0.3, 0.1); this is for drawText
    //QVector3D nearEnd(0,0,-1), farEnd(0,0,1);
    nearEnd.setX(loc1.x() + (-1-loc1.z())*(loc0.x()-loc1.x())/(loc0.z()-loc1.z()));
    nearEnd.setY(loc1.y() + (-1-loc1.z())*(loc0.y()-loc1.y())/(loc0.z()-loc1.z()));
    nearEnd.setZ(-1);
    farEnd.setX(loc1.x() + (1-loc1.z())*(loc0.x()-loc1.x())/(loc0.z()-loc1.z()));
    farEnd.setY(loc1.y() + (1-loc1.z())*(loc0.y()-loc1.y())/(loc0.z()-loc1.z()));
    farEnd.setZ(1);

    nearEnd = view_pos_to_volume_pixel_pos(nearEnd);
    farEnd = view_pos_to_volume_pixel_pos(farEnd);
}
/**
 * @brief process_right_button_hit: pop-up menu to process cell traces， MainWindow *_mainwindow
 */
int RayCastCanvas::process_right_button_hit(QMouseEvent *event){
    if (!bShowTrackResult) return 0;
    if (bWait4RightClickOnCell){ // waiting for user's annotation
        QVector3D farend, nearend;
        userClick2volumeLoc(event, nearend, farend);
        size_t user_click_cell_center_idx = cellTracker->endpoint2loc(nearend, farend, *cellSegmenter,
                                                         curr_timePoint_in_canvas);
        if(canvasRightClickOperation == RRIGHT_CLICK_OPR_EXTEND_EXIST_TRACE && curr_trace_id_from_click>=0){
            cellTracker->extendTraceWithOneAnnotation(curr_trace_id_from_click, user_click_cell_center_idx,
                                                      *cellSegmenter, curr_timePoint_in_canvas,
                                                      data_importer->filelist.at(curr_timePoint_in_canvas));
        }else if(canvasRightClickOperation == RRIGHT_CLICK_OPR_SELECT_ONE_TRACE){
            curr_trace_id_from_click = cellTracker->loc2traceId(user_click_cell_center_idx, *cellSegmenter, curr_timePoint_in_canvas,
                                                            data_importer->filelist.at(curr_timePoint_in_canvas));
            this->import_traces(curr_timePoint_in_canvas);
        }else{
            qDebug() << "Non-defined operation.";
        }
        bWait4RightClickOnCell = false; // annotation obtained
        update();
    }else{
        QList<QAction*> listAct;
        QAction *act=0, *actShowSingleCellTrace=0, *actContinueOneTrace;

        listAct.append(act = new QAction("", this->parent())); act->setSeparator(true);
        listAct.append(actShowSingleCellTrace = new QAction("Right-click to select a cell to track", this->parent()));
        listAct.append(actContinueOneTrace = new QAction("Right-click to continue the cell trace", this->parent()));

        QMenu menu;
        foreach (QAction* a, listAct) {  menu.addAction(a); }
        //menu.setWindowOpacity(POPMENU_OPACITY); // no effect on MAC? on Windows cause blink
        act = menu.exec(QCursor::pos());
        if (act==0) 	return 0;
        else if (act == actShowSingleCellTrace)
        {
            bWait4RightClickOnCell = true;
            canvasRightClickOperation = RRIGHT_CLICK_OPR_SELECT_ONE_TRACE;
            qDebug() << "actShowSingleCellTrace";
        }
        else if (act == actContinueOneTrace)
        {
            bWait4RightClickOnCell = true;
            canvasRightClickOperation = RRIGHT_CLICK_OPR_EXTEND_EXIST_TRACE;
            //bExtendExistingTrace = true;
            qDebug() << "actContinueOneTrace";
        }else{
            qDebug() << "No valid selection";
            return 0;
        }
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
void RayCastCanvas::setVolumeTimePoint(int t_at_curr_loaded_data)//,
{
    //qDebug("V3dR_GLWidget::setVolumeTimePoint = %d", t);
    if (t_at_curr_loaded_data<0) t_at_curr_loaded_data = 0;
    if (data_importer && t_at_curr_loaded_data>=data_importer->image4d->getTDim()){
        t_at_curr_loaded_data = data_importer->image4d->getTDim()-1;
    }
    this->setVolume(t_at_curr_loaded_data);

    /*! t_at_curr_loaded_data may not equals to current time point. We may only load part of the data */
    curr_timePoint_in_canvas = data_importer->curr_start_file_id + t_at_curr_loaded_data; // default is 0
    //emit changeVolumeTimePoint(t); //need?
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
    update();
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
        QVector3D originPt = m_modelViewProjectionMatrix * QVector3D(-1, 1, -1);
        // x-axis
        QVector3D newEndPt = m_modelViewProjectionMatrix * QVector3D(exceed_extent, 1, -1);
        QColor c = QColor(255,0,0);
        this->drawLine(&painter, c, QPointF(originPt.x(), originPt.y()),
                       QPointF(newEndPt.x(), newEndPt.y()), 2);
        this->drawText(&painter, c, QPointF(newEndPt.x(), newEndPt.y()), "X");
        // y-axis
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, -exceed_extent, -1);
        c = QColor(0,255,0);
        this->drawLine(&painter, c, QPointF(originPt.x(), originPt.y()),
                       QPointF(newEndPt.x(), newEndPt.y()), 2);
        this->drawText(&painter, c, QPointF(newEndPt.x(), newEndPt.y()), "Y");
        // z-axis
        newEndPt = m_modelViewProjectionMatrix * QVector3D(-1, 1, exceed_extent);
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


void RayCastCanvas::draw_makers(){
    QPainter painter(this);
    for(auto &marker : markers){
        if(marker.time_point == curr_timePoint_in_canvas){
            QVector3D Z(0, 0, 0); // instead of 0 for x and y i need worldPosition.x() and worldPosition.y() ....
            Z = Z.project(marker.ViewMatrix*marker.ModelMatrix,
                          marker.ProjectionMatrix, QRect(0, 0, marker.canvas_width, marker.canvas_height));
            QVector3D loc0 =
                    QVector3D(marker.canvas_pos.x(), marker.canvas_height - marker.canvas_pos.y(),
                    Z.z()).unproject(marker.ViewMatrix*marker.ModelMatrix, marker.ProjectionMatrix, QRect(0, 0, marker.canvas_width, marker.canvas_height));
            qDebug() << loc0;

            Z = QVector3D(0, 0, 1);
            Z = Z.project(marker.ViewMatrix*marker.ModelMatrix,
                          marker.ProjectionMatrix, QRect(0, 0, marker.canvas_width, marker.canvas_height));
            QVector3D loc1 =
                    QVector3D(marker.canvas_pos.x(), marker.canvas_height - marker.canvas_pos.y(),
                    Z.z()).unproject(marker.ViewMatrix*marker.ModelMatrix, marker.ProjectionMatrix, QRect(0, 0, marker.canvas_width, marker.canvas_height));

            qDebug() << loc1;
            // way 1: get the exact coordinate of worldPosition, e.g.: (300,200,10)
            //QVector3D curPt = worldPosition.project(m_viewMatrix*m_modelMatrix, m_projectionMatrix, QRect(0, 0, width(), height()));
            // way 2: get the normalized location of worldPosition, e.g.: (0.5, 0.3, 0.1); this is for drawText
            QVector3D nearEnd(0,0,-1), farEnd(0,0,1);
            nearEnd.setX(loc1.x() + (-1-loc1.z())*(loc0.x()-loc1.x())/(loc0.z()-loc1.z()));
            nearEnd.setY(loc1.y() + (-1-loc1.z())*(loc0.y()-loc1.y())/(loc0.z()-loc1.z()));
            farEnd.setX(loc1.x() + (1-loc1.z())*(loc0.x()-loc1.x())/(loc0.z()-loc1.z()));
            farEnd.setY(loc1.y() + (1-loc1.z())*(loc0.y()-loc1.y())/(loc0.z()-loc1.z()));

            QColor c = QColor(255,0,0);
            QVector3D p_start = m_modelViewProjectionMatrix * nearEnd;
            this->drawText(&painter, c, QPointF(p_start.x(), p_start.y()), "back");

            c = QColor(255,255,0);
            QVector3D p_end = m_modelViewProjectionMatrix * farEnd;
            this->drawText(&painter, c, QPointF(p_end.x(), p_end.y()), "front");

            painter.setPen(QPen(c, 1));
            c = QColor(0,255,0);
            this->drawLine(&painter, c, QPointF(p_start.x(), p_start.y()), QPointF(p_end.x(), p_end.y()));
            //painter.setPen(QPen(Qt::white, 3));
            //this->drawLine(&painter, Qt::white, QPointF(0, 0), QPointF(1, 1));

            //this->drawText(&painter, c, QPointF(p_end.x(), p_end.y()), "front");
        }
    }
    painter.end();
}

/**
 * @brief import_traces
 * @param movieInfo: celltracking results
 * @param t: neglect the traces that are un-related with frame t; Note this t here is not influenced by
 * the loaded data. Even we only load part of the data, as long as the trace infor in intact, t is for
 * the whole data.
 */
void RayCastCanvas::import_traces(int t){
    traces.clear();
    traces.resize(cellTracker->movieInfo.tracks.size());
    curr_timePoint_in_canvas = t;
    if(curr_trace_id_from_click < 0){ // display all valid traces in the FOV
        for(int j=0; j<cellTracker->movieInfo.tracks.size(); j++){
            if(cellTracker->movieInfo.tracks[j].size()<=5){ // if the trace has stopped before time t
                continue;
            }
            int end_time = cellTracker->movieInfo.frames[*cellTracker->movieInfo.tracks[j].rbegin()];
            int start_time = cellTracker->movieInfo.frames[*cellTracker->movieInfo.tracks[j].begin()];
            if(end_time < t || start_time > t){ // if the trace has stopped before time t
                continue;
            }

            for(auto idx : cellTracker->movieInfo.tracks[j]){
                if(cellTracker->movieInfo.frames[idx] > t) break;
                traces[j].emplace_back(QVector3D(cellTracker->movieInfo.xCoord[idx], cellTracker->movieInfo.yCoord[idx],
                                                 cellTracker->movieInfo.zCoord[idx]));
            }
        }
    }else{ // display one trace only
        for(auto idx : cellTracker->movieInfo.tracks[curr_trace_id_from_click]){
            if(cellTracker->movieInfo.frames[idx] > t) break;
            traces[curr_trace_id_from_click].emplace_back(QVector3D(cellTracker->movieInfo.xCoord[idx],
                                                           cellTracker->movieInfo.yCoord[idx],
                                                           cellTracker->movieInfo.zCoord[idx]));
        }
    }
    //traces.resize(trace_cnt);
}
/**
 * @brief draw_traces: draw cell traces
 */
void RayCastCanvas::draw_traces(){
    QPainter painter(this);
    int line_width = 3;
    for(int tr = 0; tr<traces.size(); tr++){
        auto &trace = traces[tr];
        if(trace.empty()) continue;
        cv::Vec3b cur_cl = colormap4tracking_res->at<cv::Vec3b>(tr);
        QColor c(cur_cl(0), cur_cl(1), cur_cl(2));
        //QColor c = QColor(255,0,0);
        //painter.setPen(QPen(c, 3));
        //        qDebug() << c;
        if(trace.size() > 1){
            QVector3D trace_head = volume_pixel_pos_to_view_pos(trace[0]);
            trace_head = m_modelViewProjectionMatrix * trace_head;
            this->drawPoint(&painter, c, QPointF(trace_head.x(), trace_head.y()), line_width+1);
        }
        for(int i=1; i<trace.size(); i++){
            QVector3D start_p = volume_pixel_pos_to_view_pos(trace[i-1]);
            QVector3D end_p = volume_pixel_pos_to_view_pos(trace[i]);
            //qDebug() << trace[i-1] << trace[i];
            //qDebug() << start_p << end_p;
            QVector3D p_start = m_modelViewProjectionMatrix * start_p;
            QVector3D p_end = m_modelViewProjectionMatrix * end_p;

            //qDebug() << p_start << p_end;
            //this->drawText(&painter, c, QPointF(p_start.x(), p_start.y()), "x");
            //this->drawText(&painter, c, QPointF(p_end.x(), p_end.y()), "o");
            this->drawLine(&painter, c, QPointF(p_start.x(), p_start.y()), QPointF(p_end.x(), p_end.y()), line_width);
            //painter.drawLine(QPointF(p_start.x(), p_start.y()), QPointF(p_end.x(), p_end.y()));            //flag = true;
            //break;
            this->drawPoint(&painter, c, QPointF(p_end.x(), p_end.y()), line_width+1);
        }
        //if(flag) break;
    }
    painter.end();
}
