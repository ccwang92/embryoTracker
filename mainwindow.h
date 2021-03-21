#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QScrollArea>
#include <QScrollBar>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QMenuBar>
#include <QFileDialog>
#include <QStatusBar>
//#include "emt_glwidget.h"
#include "data_importer.h"
#include "myglwidget.h"
#include "raycasting/raycastcanvas.h"
#include "cellsegmentation/cellsegment_main.h"
#include "celltracking/celltracking_main.h"

enum widget_type_choice {my_simple_test_type = 0, raycast_type = 1, vaa3d_type = 2};


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();
signals:
    void signalDataLoaded(); //data loaded sucessfully

public slots:
    void importImageSeries();
    void initControlWidgetValues();
    void updateControlPanel();
    void sendData4Segment();
    void sendData4Track();
    void transferRGBAVolume(int t);
    // parameters
public:
    QString win_title;
    QMenuBar *menuBar;
    // file menuread the image stack in
    QMenu *fileMenu;
    QAction *importImageSeriesAct;
    QAction *exitAct;
    DataImporter *data4test = 0; // functions to import data
    // edit menu
    QMenu *editMenu;
    QAction *resetViewPoint;
    QAction *bndAxesShow;
    // process menu
    QMenu *processMenu;
    QAction *segmentCell3d;
    QAction *trackCell3d;
    // control widgets on the right size
    //QGridLayout *rightSideControlLayout;
    QScrollBar *contrastScrollBar, *thresholdScrollBar;
    // major widget
    QScrollArea *glWidgetArea = 0;
    int widget_type = raycast_type;
    QGroupBox* grpBox4display_canvas;
    //EmT_GLWidget *glWidget = 0;
    MyGLWidget *glWidget_simple = 0;
    RayCastCanvas *glWidget_raycast = 0;
    // time (frame) control
    QScrollBar *timeSlider = 0;
    //QAction *exitAct;
    void createControlWidgets();
    void connectSignal();
    QScrollBar *createContrastSlider();
    QScrollBar *createThreshodSlider();
public:
    bool algorithmDebug = false;
    //QString debugDataPath = QString("/home/ccw/Desktop/embryo_res_folder/crop_embryo_data_500x500x30x40/images_downsample/1.tif");//_5f
    QString debugDataPath = QString("/home/ccw/Desktop/embryo_res_folder/"
                                    "downsample_crop_embryo_data_470x350x250x50/embryo_TM481.tif");
    //QString debugDataPath = QString("/home/ccw/Desktop/test_ims/cropped_16bit4speed/embryo_TM481.tif");
    QAction * debugButton;
    cellSegmentMain *cellSegmenter = nullptr;
    cellTrackingMain *cellTracker = nullptr;
    bool tracking_result_exist = false;
    Mat3b colormap4tracking_res;
public slots:
    void debugAlgorithm();


};
#endif // MAINWINDOW_H
