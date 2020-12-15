#pragma once

#include <functional>
#include <vector>

#include <QtMath>

#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLShaderProgram>
#include <QMessageBox>
#include "mesh.h"
#include "raycastvolume.h"
#include "trackball.h"
//#include "vtkvolume.h"
#include "../data_importer.h"

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
        qDebug("   catch: rgbaBuf = @%0p", rgbaBuf); \
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
        data_importer = _data_importer;
        setVolume();
    }
    void setVolume(long frame4display = 0);
    void setThreshold(const double threshold) {
        auto range = m_raycasting_volume ? getRange() : std::pair<double, double>{0.0, 1.0};
        m_threshold = threshold / (range.second - range.first);
        update();
    }

    void setMode(const QString& mode) {
        m_active_mode = mode;
        update();
    }

    void setBackground(const QColor& colour) {
        m_background = colour;
        update();
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

    void handleKeyPressEvent(QKeyEvent * event); //for hook to MainWindow
    void handleKeyReleaseEvent(QKeyEvent * event); //for hook to MainWindow

signals:
    void changeVolumeTimePoint(int);
public slots:
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * event);
    virtual void setVolumeTimePoint(int t);
    virtual void setLightPositionZero();
    virtual void setContrast(int relative_contrast/*[-100:100]*/);
    virtual void keyPressEvent(QKeyEvent *e){handleKeyPressEvent(e);}
    virtual void keyReleaseEvent(QKeyEvent *e){handleKeyReleaseEvent(e);}
    virtual void setBnfAxesOnOff();
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

public:
    float *depth_buffer;
    RGBA8 *total_rgbaBuf, *rgbaBuf;  // this will be updated when needs rendering (not quite sure why Vaa3d needs two vectors)
    //float sampleScale[5];
    V3DLONG start1, start2, start3, start4, start5;
    V3DLONG size1, size2, size3, size4, size5;
    V3DLONG dim1, dim2, dim3, dim4, dim5;
    V3DLONG bufSize[5]; //(x,y,z,c,t) 090731: add time dim

private:
    DataImporter *data_importer;
    QMatrix4x4 m_viewMatrix;
    QMatrix4x4 m_modelViewProjectionMatrix;
    QMatrix3x3 m_normalMatrix;
    // m_fov is the maximum vertical angle of cammera, it defines how large the fov will be
    const GLfloat m_fov = 30.0f;                                          /*!< Vertical field of view. */
    const GLfloat m_focalLength = 1.0 / qTan(M_PI / 180.0 * m_fov / 2.0); /*!< Focal length. */
    GLfloat m_aspectRatio;                                                /*!< width / height */
    GLboolean consider_transparency = false;
    QVector2D m_viewportSize;
    QVector3D m_rayOrigin; /*!< Camera position in model space coordinates. */

    QVector3D m_lightPosition {3.0, 0.0, 3.0};    /*!< In camera coordinates. */
    QVector3D m_diffuseMaterial {1.0, 1.0, 1.0};  /*!< Material colour. */
    GLfloat m_stepLength;                         /*!< Step length for ray march. */
    // no use
    GLfloat m_threshold = 50;                     /*!< Isosurface intensity threshold. */

    QColor m_background;                          /*!< Viewport background colour. */

    GLfloat m_gamma = 1.0f;//2.2 /*!< Gamma correction parameter. */

    RayCastVolume *m_raycasting_volume;
    //QPainter *canvas_painter;
    std::map<QString, QOpenGLShaderProgram*> m_shaders;
    std::map<QString, std::function<void(void)>> m_modes;
    QString m_active_mode;

    TrackBall m_trackBall {};       /*!< Trackball holding the model rotation. */
    TrackBall m_scene_trackBall {}; /*!< Trackball holding the scene rotation. */

    //center shift
    QPointF *centerShift = new QPointF(0.0, 0.0);
    GLint m_distExp = -200;

    GLuint scaled_width();
    GLuint scaled_height();

    void raycasting(const QString& shader);

    QPointF pixel_pos_to_view_pos(const QPointF& p);
    void create_noise(void);
    void add_shader(const QString& name, const QString& vector, const QString& fragment);
public:
    // rendering text
    void renderText(double x, double y, double z, QString text);
    inline GLint project(GLdouble objx, GLdouble objy, GLdouble objz,
                        const GLdouble model[16], const GLdouble proj[16],
                        const GLint viewport[4],
                        GLdouble * winx, GLdouble * winy, GLdouble * winz);
    inline void transformPoint(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]);

    // add bounding box and x-, y-, z-axes
    bool bShowAxes = true, bShowBoundingBox = false;
    BoundingBox* posXTranslateBB=0;
    BoundingBox* negXTranslateBB=0;
    BoundingBox* posYTranslateBB=0;
    BoundingBox* negYTranslateBB=0;
    BoundingBox boundingBox = 0;
    int bShowXYTranslateArrows = 0,iPosXTranslateArrowEnabled = 0,
    iNegXTranslateArrowEnabled = 0,iPosYTranslateArrowEnabled = 0,
    iNegYTranslateArrowEnabled = 0;
    RGBA32f color_line = XYZW(.5f,.5f,.7f, 1);


    void draw_bbox();
    void drawString(float x, float y, float z, const char* text, int shadow=0, int fontsize=0);
    void setBoundingBoxSpace(BoundingBox BB);
    virtual void getBoundingBox(BoundingBox& bb) {bb = boundingBox;};
    virtual void drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth=1, float AlineWidth=3);
};
