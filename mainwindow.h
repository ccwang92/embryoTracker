#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//#include "emt_glwidget.h"
#include "data_importer.h"
#include "myglwidget.h"
#include "raycasting_renderer/src/raycastcanvas.h"
#include <QMainWindow>
#include <QScrollArea>
#include <QScrollBar>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QMenuBar>
#include <QFileDialog>
#include <QStatusBar>

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
    void updateControlPanel();


    // parameters
public:
    QString win_title;
    // read the image stack in
    QMenuBar *menuBar;
    QMenu *fileMenu;
    QAction *importImageSeriesAct;
    QAction *exitAct;
    DataImporter *data4test = 0; // functions to import data

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
    //void

};
#endif // MAINWINDOW_H
