#include "mesh.h"


/*!
 * \brief Create a vertex array object (VAO) with given vertices and indices.
 * \param gl_context OpenGL context for the buffers.
 * \param vertices Vector of vertices, as strided x-y-z coordinates.
 * \param indices Indices for the mesh.
 */
Mesh::Mesh(const std::vector<GLfloat>& vertices, const std::vector<GLuint> &indices)
    : m_vertices_count {vertices.size()}
    , m_indices_count {indices.size()}
{
    /** QT style ***/
    this->initializeOpenGLFunctions();
    m_vao = new QOpenGLVertexArrayObject();
    m_vertex_VBO = new QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
    m_index_VBO = new QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);

    m_vao->create();
    m_vertex_VBO->create();
    m_index_VBO->create();

    m_vao->bind();
        m_vertex_VBO->bind();
        m_vertex_VBO->setUsagePattern(QOpenGLBuffer::StaticDraw);
        m_vertex_VBO->allocate(vertices.data(), vertices.size() * sizeof(vertices[0]));


        m_index_VBO->bind();
        m_index_VBO->setUsagePattern(QOpenGLBuffer::StaticDraw);
        m_index_VBO->allocate(indices.data(),indices.size()*sizeof(indices[0]));


        //Vertex position
        glEnableVertexAttribArray(POSITION);
        glVertexAttribPointer(POSITION, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
//        program->enableAttributeArray(POSITION);
//        program->setAttributeBuffer(0,GL_FLOAT,0,sizeof(Vertex));

//        program->enableAttributeArray(1);
//        program->setAttributeBuffer(1,GL_FLOAT,offsetof(Vertex,Normal),3,sizeof(Vertex));

//        glEnableVertexAttribArray(2);
//        program->setAttributeBuffer(2,GL_FLOAT,offsetof(Vertex,TextCoords),2,sizeof(Vertex));
    m_vao->release();
    /** OPENGL style ***/
//    initializeOpenGLFunctions();
//    // Generates and populates a VBO for the vertices
//    glDeleteBuffers(1, &(m_vertex_VBO));
//    glGenBuffers(1, &(m_vertex_VBO));
//    glBindBuffer(GL_ARRAY_BUFFER, m_vertex_VBO);
//    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertices[0]), vertices.data(), GL_STATIC_DRAW);

//    // Generates and populates a VBO for the element indices
//    glDeleteBuffers(1, &(m_index_VBO));
//    glGenBuffers(1, &(m_index_VBO));
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_VBO);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(indices[0]), indices.data(), GL_STATIC_DRAW);

//    // Creates a vertex array object (VAO) for drawing the mesh
//    glDeleteVertexArrays(1, &(m_vao));
//    glGenVertexArrays(1, &(m_vao));
//    glBindVertexArray(m_vao);
//    glBindBuffer(GL_ARRAY_BUFFER, m_vertex_VBO);
//    glEnableVertexAttribArray(POSITION);
//    glVertexAttribPointer(POSITION, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_index_VBO);
//    glBindVertexArray(0);
}


/*!
 * \brief Destructor.
 */
Mesh::~Mesh()
{
}


/*!
 * \brief Render the mesh.
 */
void Mesh::paint(void)
{
    //qDebug("is vao created ? %d", )
    m_vao->bind();
        //glBindVertexArray(m_vao);
        glDrawElements(GL_TRIANGLES, m_indices_count, GL_UNSIGNED_INT, nullptr);
        //glBindVertexArray(0);
    m_vao->release();
}
