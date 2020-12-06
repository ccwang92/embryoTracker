#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "emt_glwidget.h"
#include "data_importer.h"
#include <QScrollArea>
#include <QScrollBar>

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow();
    ~MainWindow();
public slots:
    void importImageSeries();
//public signals:


    // parameters
public:
    QString win_title;

    // major widget
    //EmT_GLWidget glWidget;

    // time (frame) control
    QScrollArea *glWidgetArea;
    QScrollBar *timeSlider;

    // read the image stack in
    QMenu *fileMenu;
    QAction *importImageSeriesAct;
    QAction *exitAct;

    DataImporter *data4test = 0;
};
#endif // MAINWINDOW_H
