#ifndef EMT_GLWIDGET_H
#define EMT_GLWIDGET_H

//#include "src_3rd/3drenderer/GLee2glew.h"
#include "data_importer.h"
//#include "renderer_emt.h"

#include <QOpenGLWidget>
#include <QOpenGLFunctions>



#define ANGLE_X0 (15)			//degree
#define ANGLE_Y0 (360-20)		//degree
#define ANGLE_Z0 (360-2)		//degree
#define ANGLE_TICK 1 			//divided
#define MOUSE_SENSITIVE 1.0f
#define SHIFT_RANGE 100 		//percent of bound
#define ZOOM_RANGE  100         //percent of fov
#define ZOOM_RANGE_RATE 5       //zoom rate of fov
#define CLIP_RANGE  200 		//size of (-100,100)
#define ZTHICK_RANGE 20			//times
#define TRANSPARENT_RANGE 100   //nonlinear divided

#define POPMENU_OPACITY 1

// 081025 by RZC
#define WIN_SIZEX 1024 //800
#define WIN_SIZEY 768  //800
#define CTRL_SIZEX 350
#define MINVIEW_SIZEX 700  //800
#define MINVIEW_SIZEY 700  //800

#define QEvent_Ready (QEvent::User +1)
#define QEvent_OpenFiles (QEvent::User +2)
#define QEvent_DropFiles (QEvent::User +3)
#define QEvent_InitControlValue (QEvent::User +4)
#define QEvent_OpenLocal3DView (QEvent::User +5)
#define QEvent_HistoryChanged (QEvent::User +6) //20170801 RZC: a notify after Undo Histoty changed


class EmT_GLWidget : public QOpenGLWidget
{
public:
    EmT_GLWidget(DataImporter* idep, QWidget* mainWindow=0, QString title="");
    ~EmT_GLWidget();

    virtual void deleteRenderer();  //090710 RZC: to delete renderer before ~V3dR_GLWidget()
    virtual void createRenderer();  //090710 RZC: to create renderer at any time

    DataImporter* getiDrawExternalParameter() {return _idep;}
    QWidget* getMainWindow() {return mainwindow;}
    //Renderer* getRenderer()   {return renderer;}
    //const Renderer* getRenderer() const {return renderer;} // const version CMB
    QString getDataTitle()    {return data_title;}
    void setDataTitle(QString newdt) {data_title = newdt;}

protected:
    virtual void initializeGL();
    virtual void resizeGL(int width, int height);
    virtual void paintGL();

//    virtual void paintEvent(QPaintEvent *event);
//    virtual void focusInEvent(QFocusEvent* e);
//    virtual void focusOutEvent(QFocusEvent* e);
//    virtual void enterEvent(QEvent *e);
//    virtual void leaveEvent(QEvent *e);
//    virtual void mousePressEvent(QMouseEvent *event);
//    virtual void mouseReleaseEvent(QMouseEvent *event);
//    virtual void mouseMoveEvent(QMouseEvent *event);
//    virtual void wheelEvent(QWheelEvent *event);
//    //virtual void mouseDoubleClickEvent ( QMouseEvent * event ) {};

//    //virtual void keyPressEvent(QKeyEvent * e) {handleKeyPressEvent(e);}
//    //virtual void keyReleaseEvent(QKeyEvent * e) {handleKeyReleaseEvent(e);}

//    virtual void closeEvent(QCloseEvent* e); //for self closing
//    virtual bool event(QEvent* e);       //090427 RZC:  for QHelpEvent of ToolTip
//    virtual void customEvent(QEvent* e); // for QEvent_OpenFiles, by RZC 081002

    DataImporter* _idep;
    QWidget *mainwindow;
    Renderer_EmT* renderer;
    QString data_title;

public:
    int _renderMode;
    //unsigned char * data;
    int _data_size[5];

    char tipBuf[1000];
    bool _holding_num[10];

    int viewW, viewH;
    GLdouble mRot[16];
    GLdouble mAltC[4]; //0806 RZC: alternate rotation center (for non-center rotation)
    int alt_rotation;
    static const int flip_X= +1, flip_Y= -1, flip_Z= -1; // make y-axis downward conformed with image coordinate
    QPoint lastPos;

    float _xRot, _yRot, _zRot, dxRot, dyRot, dzRot;
    // int _zoom, _xShift, _yShift, _zShift, dxShift, dyShift, dzShift;
    float _zoom, _xShift, _yShift, _zShift, dxShift, dyShift, dzShift; // CMB 2011 Feb 07
    int _xCut0, _xCut1, _yCut0, _yCut1, _zCut0, _zCut1, _fCut;
    int dxCut, dyCut, dzCut, lockX, lockY, lockZ;
    int _xCS, _yCS, _zCS;
    int _xClip0, _xClip1, _yClip0, _yClip1, _zClip0, _zClip1;
    int _CStransparency, _markerSize, _curChannel;
    float _thickness;
    int _Bright, _Contrast, sUpdate_bright, sUpdate_track;
    bool _showAxes, _showBoundingBox, _absRot, _orthoView, _clipBoxEnable;
    bool _volCompress, _volFilter;

    //RGBA8 backgroundColor; // record current non-black backgroundColor

    int _volumeTimePoint; float volumeTimPoint_fraction;
    void init_members()
    {
        _renderMode = 0;
        for (int i=0; i<5; i++)	_data_size[i] = 0;

        for (int i=0; i<10; i++) _holding_num[i] = false;

        viewW=viewH=0;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                mRot[i*4 +j] = ((i==j)? 1 : 0); // Identity matrix
        for (int i=0; i<4; i++)
                mAltC[i] = ((i==3)? 1 : 0);		// Original point
        alt_rotation =0;

        _xRot=_yRot=_zRot= dxRot=dyRot=dzRot=
        _zoom=_xShift=_yShift=_zShift= dxShift=dyShift=dzShift=
        _xCut0=_xCut1=_yCut0=_yCut1=_zCut0=_zCut1=_fCut=
        dxCut=dyCut=dzCut= lockX=lockY=lockZ=
        _xCS=_yCS=_zCS=
        _xClip0=_xClip1=_yClip0=_yClip1=_zClip0=_zClip1 =0;
        _thickness =1;
        _CStransparency=0; _markerSize=1, _curChannel=-1;

        _Bright=_Contrast=sUpdate_bright=sUpdate_track=0;
        _showAxes = _showBoundingBox = _absRot = _orthoView =false;
        _clipBoxEnable = _volCompress = _volFilter =true;

        _volumeTimePoint=0;
        volumeTimPoint_fraction=0;
    }
protected:
    scene::GLscene *proxty;
//public slots:
//signals:

};

#endif // EMT_GLWIDGET_H
