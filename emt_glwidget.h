#ifndef EMT_GLWIDGET_H
#define EMT_GLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>

class EmT_GLWidget : public QOpenGLWidget
{
public:
    //EmT_GLWidget(QWidget *parent) : QOpenGLWidget(parent) { };
    EmT_GLWidget();
protected:
    void initializeGL() override
    {
        // Set up the rendering context, load shaders and other resources, etc.:
        QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
        f->glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    }

    void resizeGL(int w, int h) override
    {
        // Update projection matrix and other size related settings:
        //m_projection.setToIdentity();
        //m_projection.perspective(45.0f, w / float(h), 0.01f, 100.0f);
    }

    void paintGL() override
    {
        // Draw the scene:
        QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();
        f->glClear(GL_COLOR_BUFFER_BIT);

    }

};

#endif // EMT_GLWIDGET_H
