#pragma once

#include <functional>
#include <vector>

#include <QtMath>
#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QMessageBox>
#include <math.h>
#include "mesh.h"
#include "raycastvolume.h"
#include "trackball.h"
#include <opencv2/core.hpp>  //mat4b
#include "../data_importer.h"
#include "../cellsegmentation/types_define.h" // MarkerPos
#include "../cellsegmentation/cellsegment_main.h"
#include "../celltracking/celltracking_main.h"
// if error then close
// clean memory before MessageBox, otherwise MessageBox maybe could not be created correctly
#define ERROR_MessageBox(title, type, what) { \
    cleanData(); \
    QMessageBox::critical( 0, title, QObject::tr("%1: OUT OF MEMORY or GL ERROR.\n---%2 exception: %3")\
            .arg(title).arg(type).arg(what) + "\n\n" + \
        QObject::tr("Please close some images or views to release memory, then try again.\n\n") ); \
}
#define CATCH_handler( func_name ) \
    catch (std::exception & e) { \
        \
        ERROR_MessageBox( func_name, "std", e.what() ); \
        return; \
        \
    } catch (const char* s) { \
        \
        ERROR_MessageBox( func_name, "GL", s ); \
        return; \
        \
    } catch (...) { \
        \
        ERROR_MessageBox( func_name, "UNKOWN", "unknown exception" ); \
        return; \
        \
    }

#define DELETE_AND_ZERO(p)	{ if ((p)!=NULL) delete (p); (p) = NULL; }
/*!
 * \brief Class for a raycasting canvas widget.
 */
class RayCastCanvas : public QOpenGLWidget, protected QOpenGLExtraFunctions
{
    Q_OBJECT
public:
    explicit RayCastCanvas(QWidget *parent = nullptr);
    ~RayCastCanvas();
    virtual void cleanData();
    void setStepLength(const GLfloat step_length) {
        m_stepLength = step_length;
        update();
    }
    void initVolume(DataImporter *_data_importer) {
        //m_raycasting_volume = new RayCastVolume();
        //m_raycasting_volume->initMesh();
        //m_raycasting_volume->create_noise();
        markers.clear();
        data_importer = _data_importer;

        setVolume();
    }
    void setVolume(long frame4display = 0);
    void setVolumeWithMask(long frame4display, unsigned char* mask);
    void setThreshold(const double threshold) {
        auto range = m_raycasting_volume ? getRange() : std::pair<double, double>{0.0, 1.0};
        m_threshold = threshold / (range.second - range.first);
        update();
    }

    void setMode(const QString& mode) {
        m_active_mode = mode;
        update();
    }
    void resetMode(void){
        m_active_mode = m_init_mode;
        update();
    }
    void setBackground(const QColor& colour) {
        m_background = colour;
        update();
    }
    QString getMode(void) {
        return m_active_mode;
    }
    std::vector<QString> getModes(void) {
        std::vector<QString> modes;
        for (const auto& [key, val] : m_modes) {
            modes.push_back(key);
        }
        return modes;
    }

    QColor getBackground(void) {
        return m_background;
    }

    std::pair<double, double> getRange(void) {
        return m_raycasting_volume->range();
    }

    DataImporter* getDataImporter(){return data_importer;}
    RayCastVolume* getRenderer() {return m_raycasting_volume;}
    void handleKeyPressEvent(QKeyEvent * event); //for hook to MainWindow
    void handleKeyReleaseEvent(QKeyEvent * event); //for hook to MainWindow

signals:
    void changeVolumeTimePoint(int);
    //void changeBnfAxesOnOff(bool);
public slots:
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * event);
    virtual void setVolumeTimePoint(int t);
    virtual void setLightPositionZero();
    virtual void setContrast(int relative_contrast/*[-100:100]*/);
    virtual void setThreshold(int intensity_threshold);
    virtual void setRangeXMIN(int xmin);
    virtual void setRangeYMIN(int ymin);
    virtual void setRangeZMIN(int zmin);
    virtual void setRangeXMAX(int xmax);
    virtual void setRangeYMAX(int ymax);
    virtual void setRangeZMAX(int zmax);
    virtual void keyPressEvent(QKeyEvent *e){handleKeyPressEvent(e);}
    virtual void keyReleaseEvent(QKeyEvent *e){handleKeyReleaseEvent(e);}
    virtual void setBnfAxesOnOff();
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    //void paintEvent(QPaintEvent *event);

public:
    float *depth_buffer;
    bool flag_rgba_display = false;
    RGBA8 *total_rgbaBuf {0}, *rgbaBuf{0};  // this will be updated when needs rendering (not quite sure why Vaa3d needs two vectors)
    //float sampleScale[5];
    //V3DLONG start1, start2, start3, start4, start5;
    //V3DLONG size1, size2, size3, size4, size5;
    //V3DLONG dim1, dim2, dim3, dim4, dim5;
    V3DLONG bufSize[5]; //(x,y,z,c,t) 090731: add time dim
    int curr_timePoint_in_canvas = -1;

public:
    // visualize the tracking results
    //MainWindow *mainwindow {nullptr};
    bool bShowTrackResult = false;
    bool bShowSingleTrace = false, bWait4RightClickOnCell = false;
    int curr_trace_id_from_click = -1;
    cv::Mat4b rgb_frame = cv::Mat(); // colorful data with original data overlaid with traces
    std::vector<std::vector<QVector3D>> traces;
    cv::Mat3b *colormap4tracking_res {0};
    cellSegmentMain *cellSegmenter = nullptr;
    cellTrackingMain *cellTracker = nullptr;
    void draw_traces();
    void import_traces(int t); // get traces
    QVector3D userClick2volumeLoc(QMouseEvent *event);
    int process_right_button_hit(QMouseEvent *event); //MainWindow *_mainwindow
public:
    // visualize the markers
    bool bShowMarkers = true;
    std::vector<MarkerPos> markers;
    // functions for visualize markers
    // given a marker, find the starting and ending point of the corresponding ray intersecting the data
    void MarkerPos_to_NearFarPoint(const MarkerPos & marker, QVector3D &loc0, QVector3D &loc1);
    void draw_makers();
private:
    DataImporter *data_importer {0};
    QMatrix4x4 m_viewMatrix;
    QMatrix4x4 m_modelMatrix;
    QMatrix4x4 m_projectionMatrix;
    QMatrix4x4 m_modelViewProjectionMatrix;
    QMatrix3x3 m_normalMatrix;
    QVector3D m_leftup_xyz {0, 0, 0};
    QVector3D m_rightbottom_xyz {1.0, 1.0, 1.0};
    QVector2D m_viewportSize;
    QVector3D m_rayOrigin; /*!< Camera position in model space coordinates. */
    // m_fov is the maximum vertical angle of cammera, it defines how large the fov will be
    const GLfloat m_fov = 30.0f;                                          /*!< Vertical field of view. */
    const GLfloat m_focalLength = 1.0 / qTan(M_PI / 180.0 * m_fov / 2.0); /*!< Focal length. */
    GLfloat m_aspectRatio;                                                /*!< width / height */
    GLboolean m_consider_transparency = false;
    GLfloat m_min_valid_intensity = 0;


    QVector3D m_lightPosition {3.0, 0.0, 3.0};    /*!< In camera coordinates. */
    QVector3D m_diffuseMaterial {1.0, 1.0, 1.0};  /*!< Material colour. */
    GLfloat m_stepLength;                         /*!< Step length for ray march. */
    // no use
    GLfloat m_threshold = 50;                     /*!< Isosurface intensity threshold. */

    QColor m_background;                          /*!< Viewport background colour. */
    bool m_gamma_init = false;
    float m_gamma0 = 1.0f;
    GLfloat m_gamma = 1.0f;//2.2 /*!< Gamma correction parameter. */

    RayCastVolume *m_raycasting_volume;
    //QPainter *canvas_painter;
    std::map<QString, QOpenGLShaderProgram*> m_shaders;
    std::map<QString, std::function<void(void)>> m_modes;
    QString m_active_mode, m_init_mode;

    TrackBall m_trackBall {};       /*!< Trackball holding the model rotation. */
    TrackBall m_scene_trackBall {}; /*!< Trackball holding the scene rotation. */

    //center shift
    QPointF *centerShift = new QPointF(0.0, 0.0);
    GLint m_distExp = -200;

    GLuint scaled_width();
    GLuint scaled_height();

    void raycasting(const QString& shader);

    QVector3D view_pos_to_volume_pixel_pos(const QVector3D& p);
    QVector3D volume_pixel_pos_to_view_pos(const QVector3D& p);
    QPointF canvas_pixel_pos_to_view_pos(const QPointF& p);
    QPointF view_pos_to_canvas_pixel_pos(const QPointF& p);
    void create_noise(void);
    void add_shader(const QString& name, const QString& vector, const QString& fragment);
public:
    // Rendering axes in the figure
    bool bShowAxes = false;
    void drawInstructions(QPainter *painter);
    void drawLine(QPainter *painter, QColor c=QColor(50,50,50), QPointF p0 = QPointF(0, 0),
                   QPointF p1 = QPointF(1, 1), int lineWidth = 1);
    void drawText(QPainter *painter, QColor c=QColor(50,50,50), QPointF p = QPointF(0, 0),
                   QString text = QString("text"));
    void draw_axes();

    /** another way to draw text on volume */
    void renderText(double x, double y, double z, QString text);
    inline GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
                        const GLdouble model[16], const GLdouble proj[16],
                        const GLint viewport[4],
                        GLdouble * winx, GLdouble * winy, GLdouble * winz);
    inline void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);

    /** Vaa3d's method to add bounding box and x-, y-, z-axes */
    //    bool bShowBoundingBox = false;
    //    BoundingBox* posXTranslateBB=0;
    //    BoundingBox* negXTranslateBB=0;
    //    BoundingBox* posYTranslateBB=0;
    //    BoundingBox* negYTranslateBB=0;
    //    BoundingBox boundingBox = 0;
    //    int bShowXYTranslateArrows = 0,iPosXTranslateArrowEnabled = 0,
    //    iNegXTranslateArrowEnabled = 0,iPosYTranslateArrowEnabled = 0,
    //    iNegYTranslateArrowEnabled = 0;
    //    RGBA32f color_line = XYZW(.5f,.5f,.7f, 1);
    //    void drawString(float x, float y, float z, const char* text, int shadow=0, int fontsize=0);
    //    void setBoundingBoxSpace(BoundingBox BB);
    //    virtual void getBoundingBox(BoundingBox& bb) {bb = boundingBox;};
    //    virtual void drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth=1, float AlineWidth=3);
};
