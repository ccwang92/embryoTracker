
#pragma once

#include <QMatrix4x4>
#include <QOpenGLExtraFunctions>
#include <QVector3D>
//#include "../src_3rd/basic_c_fun/basic_4dimage.h""
#include "mesh.h"
//#include "raycastcanvas.h"
#include "../src_3rd/basic_c_fun/color_xyz.h"
/*!
 * \brief Class for a raycasting volume.
 */
class RayCastVolume : protected QOpenGLExtraFunctions
{
public:
    RayCastVolume();
    virtual ~RayCastVolume();

    void initMesh();
    void transfer_volume(void*data, double p_min, double p_max, long sx,
                         long sy, long sz, long sc);
    void create_noise(void);
    void paint(void);
    std::pair<double, double> range(void);

    void set_range(double p_min, double p_max){
        this->m_range = std::make_pair(p_min, p_max);
    }
    void set_spacing(QVector3D _m_spacing){
        this->m_spacing = _m_spacing;
    }
    void set_origin(QVector3D _m_origin){
        this->m_origin = _m_origin;
    }
    void set_size(QVector3D _m_size){
        this->m_size = _m_size;
    }
    QVector3D get_size(){
        return this->m_size;
    }
    /*!
     * \brief Get the extent of the volume.
     * \return A vector holding the extent of the bounding box.
     *
     * The extent is normalised such that the longest side of the bounding
     * box is equal to 1.
     */
    QVector3D extent(void) {
        auto e = m_size * m_spacing;
        return e / std::max({e.x(), e.y(), e.z()});
    }

    /*!
     * \brief Return the model matrix for the volume.
     * \param shift Shift the volume by its origin.
     * \return A matrix in homogeneous coordinates.
     *
     * The model matrix scales a two-unit side cube to the
     * extent of the volume.
     */
    QMatrix4x4 modelMatrix(bool shift = false) {
        QMatrix4x4 modelMatrix;
        if (shift) {
            modelMatrix.translate(-m_origin / scale_factor());
        }
        modelMatrix.scale(0.5f * extent());
        return modelMatrix;
    }

    /*!
     * \brief Top planes forming the AABB.
     * \param shift Shift the volume by its origin.
     * \return A vector holding the intercept of the top plane for each axis.
     */
    QVector3D top(bool shift = false) {
        auto t = extent() / 2.0;
        if (shift) {
            t -= m_origin / scale_factor();
        }
        return t;
    }

    /*!
     * \brief Bottom planes forming the AABB.
     * \param shift Shift the volume by its origin.
     * \return A vector holding the intercept of the bottom plane for each axis.
     */
    QVector3D bottom(bool shift = false) {
        auto b = -extent() / 2.0;
        if (shift) {
            b -= m_origin / scale_factor();
        }
        return b;
    }

private:
    GLuint m_volume_texture = 0;
    GLuint m_noise_texture = 0;
    Mesh m_cube_vao;
    std::pair<double, double> m_range = std::make_pair(0,255);
    QVector3D m_origin = QVector3D(0,0,0);
    QVector3D m_spacing = QVector3D(1,1,1);
    QVector3D m_size = QVector3D(500,500,300);

    float scale_factor(void);
public:
    void * glWidget;
    // add bounding box and x-, y-, z-axes
//    bool bShowAxes = true, bShowBoundingBox = false;
//    BoundingBox* posXTranslateBB=0;
//    BoundingBox* negXTranslateBB=0;
//    BoundingBox* posYTranslateBB=0;
//    BoundingBox* negYTranslateBB=0;
//    BoundingBox boundingBox = 0;
//    int bShowXYTranslateArrows = 0,iPosXTranslateArrowEnabled = 0,
//    iNegXTranslateArrowEnabled = 0,iPosYTranslateArrowEnabled = 0,
//    iNegYTranslateArrowEnabled = 0;
//    RGBA32f color_line = XYZW(.5f,.5f,.7f, 1);


//    void draw_bbox();
//    void drawString(float x, float y, float z, const char* text, int shadow=0, int fontsize=0);
//    void setBoundingBoxSpace(BoundingBox BB);
//    virtual void getBoundingBox(BoundingBox& bb) {bb = boundingBox;};
//    virtual void drawBoundingBoxAndAxes(BoundingBox BB, float BlineWidth=1, float AlineWidth=3);
};
