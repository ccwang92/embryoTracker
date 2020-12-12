#pragma once

#include <QOpenGLExtraFunctions>

/*!
 * \brief Class to represent a mesh.
 */
class Mesh : protected QOpenGLExtraFunctions
{
public:
    Mesh(const std::vector<GLfloat>& vertices, const std::vector<GLuint>& indices);
    virtual ~Mesh();

    void paint(void);

private:
    GLuint m_vao {0};
    GLuint m_vertex_VBO {0};
    GLuint m_normal_VBO {0};
    GLuint m_index_VBO {0};
    size_t m_vertices_count {0};
    size_t m_indices_count {0};

    enum AttributeArray {POSITION=0, NORMAL=1, TEXCOORD=2};
};
