#pragma once

//#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>
/*!
 * \brief Class to represent a mesh.
 */
class Mesh : protected QOpenGLExtraFunctions
{
public:
    Mesh(const std::vector<GLfloat>& vertices, const std::vector<GLuint>& indices);
    Mesh(){};
    virtual ~Mesh();

    void paint(void);

private:
    QOpenGLVertexArrayObject *m_vao {0};
    QOpenGLBuffer *m_vertex_VBO {0};
    QOpenGLBuffer *m_normal_VBO {0};
    QOpenGLBuffer *m_index_VBO {0};
    size_t m_vertices_count {0};
    size_t m_indices_count {0};

    enum AttributeArray {POSITION=0, NORMAL=1, TEXCOORD=2};
};
