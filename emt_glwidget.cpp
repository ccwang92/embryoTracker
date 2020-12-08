#include "emt_glwidget.h"
#include "src_3rd/basic_c_fun/color_xyz.h"
#include <QCoreApplication>
#include <QEvent>

EmT_GLWidget::EmT_GLWidget(DataImporter* idep, QWidget* mainWindow, QString title){
    qDebug("EmT_GLWidget::EmT_GLWidget ========================================");

    this->_idep = idep;
    this->mainwindow = mainWindow;
    this->data_title = title;
    //this->renderer = 0;

    //this->show_progress_bar = true;
    init_members();
    //makeCurrent(); //090729: this make sure created GL context

    //setFocusPolicy(Qt::WheelFocus); // accept KeyPressEvent when mouse wheel move, by RZC 080831
    setFocusPolicy(Qt::StrongFocus); // accept KeyPressEvent when mouse click, by RZC 081028
}
EmT_GLWidget::~EmT_GLWidget()
{
    deleteRenderer();//090711 RZC: maybe too late, because some version Qt destroyed GL context before here.
}

void EmT_GLWidget::deleteRenderer() {/*makeCurrent();*/DELETE_AND_ZERO(renderer);} //090710 RZC: to delete renderer before ~EmT_GLWidget()
void EmT_GLWidget::createRenderer() {/*makeCurrent();*/ deleteRenderer(); initializeGL();} //090710 RZC: to create renderer at any time

void EmT_GLWidget::initializeGL()
{
    qDebug("EmT_GLWidget::initializeGL");

    //test if OpenGl version is too old
    GLeeInit();
    //==============================================================================
    // OpenGL hardware supporting detection
    const char* glversion = (const char*)glGetString(GL_VERSION);
    if (strlen(glversion)>3 && glversion[0]<'2')
    {
        qDebug("   *** You OpenGL version (%s) is under 2.0, switch to Cross-Section type.", glversion);
    }
    renderer = new Renderer_EmT(this);
    //if (renderer) renderer->selectMode = Renderer_EmT::defaultSelectMode;

    // set renderer to initial state
    qDebug("EmT_GLWidget::preparingRenderer");
    if (renderer)
    {
        renderer->setupData(this->_idep);
        //if (renderer->hasError())	POST_CLOSE(this);
        renderer->getLimitedDataSize(_data_size); //for update slider size
    }

    if (renderer)
    {
        renderer->initialize(); //090705 RZC
        //if (renderer->hasError())	POST_CLOSE(this);
    }
    //=============================================================================

    // when initialize done, update status of control widgets
    //QCoreApplication::sendEvent(this, new QEvent(QEvent::Type(QEvent_InitControlValue) ));
    //QCoreApplication::postEvent(this, new QEvent(QEvent::Type(QEvent_OpenFiles) ));
    //QCoreApplication::postEvent(this, new QEvent(QEvent::Type(QEvent_Ready) ));
}


void EmT_GLWidget::resizeGL(int width, int height)
{
    //qDebug(" renderer->setupView( %d, %d )", width, height);
    viewW = width; viewH = height;
    if (renderer) renderer->setupView(width,height);
}

void EmT_GLWidget::paintGL()
{
    //if (renderer && renderer->hasError())  POST_CLOSE(this);

    //QTime qtime; qtime.start();

    //the following translation & rotation operations are carried out in view space
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //the order is important

    //GET current rotation pose by GL matrix stack
    glPushMatrix();
    {
        glLoadIdentity();
        //last absolute rotation pose
        glMultMatrixd(mRot);
        //current relative small rotation, always around center of model
        {
            XYZ R(dxRot, dyRot, dzRot);  					//qDebug("R= %f %f %f", R.x, R.y, R.z);
            dxRot=dyRot=dzRot=0;  // clear relative rotation step

            double angle = norm(R)/(float)ANGLE_TICK;       //qDebug("angle=%f", angle);
            if (angle)
            {
                normalize(R);          						//qDebug("R= %f %f %f", R.x, R.y, R.z);
                glRotated( angle,  R.x, R.y, R.z);
            }
        }
        //save current absolute rotation pose
        glGetDoublev(GL_MODELVIEW_MATRIX, mRot);
        for (int i=0; i<3; i++)
            mRot[i*4 +3]=mRot[3*4 +i]=0; mRot[3*4 +3]=1; // only reserve rotation, remove translation in mRot
    }
    glPopMatrix();

    //SET translation
    {
        //absolute translation
        XYZ T(_xShift, _yShift, _zShift);  				//qDebug("T= %f %f %f", T.x, T.y, T.z);
        //XYZ T(_xShift, _yShift, 0); // force zShift=0
        dxShift=dyShift=dzShift=0;  // clear relative shift step

        double s = 1.4/(float)SHIFT_RANGE;  // 1.4 ~ sqrt(2);
        T = T*s;
        glTranslated( T.x, T.y, T.z );
    }

    //SET current absolute rotation pose at alternate rotation center
    if (alt_rotation)	glTranslated( mAltC[0]*flip_X, mAltC[1]*flip_Y, mAltC[2]*flip_Z );
    glMultMatrixd(mRot);
    if (alt_rotation)	glTranslated( -mAltC[0]*flip_X, -mAltC[1]*flip_Y, -mAltC[2]*flip_Z );


    glScaled(flip_X,flip_Y,flip_Z); // make y-axis downward conformed with image coordinate

    //glScaled(1,1, _thickness); // here may be out of view clip-space, not used

    //=========================================================================
    // normalized space of [-1,+1]^3;
    if (renderer)
    {
        renderer->paint();
        renderer->blendBrighten(_Bright/100.f, (_Contrast+100)/100.f); // fast 8-bit precision
    }

    //qDebug("paint frame cost time = %g sec", qtime.elapsed()*0.001);

    //CHECK_GLError_print(); //090715,090723
}

