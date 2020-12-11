#include "mainwindow.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_message.h"

/**
 * brief: create main window
 * */
MainWindow::MainWindow()
{
    //QVBoxLayout *mainLayout = new QVBoxLayout();
    createControlWidgets();
    // initialize data4test as null
    data4test = new DataImporter(); //do we need?

    setMenuBar(menuBar);
    setCentralWidget(grpBox4display_canvas);
    //layout->addWidget(rightSideControlLayout);
    //setLayout(mainLayout);
    connectSignal();
}
void MainWindow::initControlWidgetValues(){
    //timeSlider no need
    contrastScrollBar->setValue(0);
}
/**
 * brief: Create the major layout of the main window
 * */
void MainWindow::createControlWidgets()
{
    menuBar = new QMenuBar;
    /************ file menu define *****************/
    fileMenu = new QMenu(tr("&File"), this);
    menuBar->addMenu(fileMenu);
    // import image stacks
    importImageSeriesAct = new QAction(tr("&Import image stack..."), this);
    importImageSeriesAct->setShortcut(tr("Ctrl+I"));
    // status has not been defined
    //importImageFileAct->setStatusTip(tr("Import general image series"));

    fileMenu->addAction(importImageSeriesAct);
    // separator
    fileMenu->addSeparator();
    // exit the program
    exitAct = new QAction(tr("E&xit"), this);
    exitAct->setShortcut(tr("Ctrl+Q"));
    //exitAct->setStatusTip(tr("Exit the application"));
    fileMenu->addAction(exitAct);
    /************ edit menu define *****************/
    editMenu = new QMenu(tr("&Edit"), this);
    menuBar->addMenu(editMenu);
    // reset view point
    resetViewPoint = new QAction(tr("&Reset ViewPoint"), this);
    resetViewPoint->setShortcut(tr("Ctrl+R"));
    // status has not been defined
    //importImageFileAct->setStatusTip(tr("Import general image series"));
    editMenu->addAction(resetViewPoint);
    /***************** About *************************/
    QMenu * aboutMenu= new QMenu(tr("&About"), this);
    menuBar->addMenu(aboutMenu);
    /*************** display grid *******************/
    grpBox4display_canvas = new QGroupBox();//tr("Canvas")
    // area to show volume
    glWidgetArea = new QScrollArea;
    glWidgetArea->setWidgetResizable(true);
    glWidgetArea->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    glWidgetArea->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    glWidgetArea->setSizePolicy(QSizePolicy::Expanding,QSizePolicy::Expanding); //QSizePolicy::Ignored, QSizePolicy::Ignored);
    glWidgetArea->setMinimumSize(500,500);//(MINVIEW_SIZEX, MINVIEW_SIZEY);
    if (widget_type == my_simple_test_type){
        glWidget_simple = new MyGLWidget(this);
        glWidgetArea->setWidget(glWidget_simple);
    }
    else if (widget_type == raycast_type){
        glWidget_raycast = new RayCastCanvas(this);
        glWidgetArea->setWidget(glWidget_raycast);
    }
    else //vaa3d_type
    {
        //
    }
    // time slider
    timeSlider = new QScrollBar(Qt::Horizontal);
    timeSlider->setRange(0,0);
    timeSlider->setSingleStep(1);
    timeSlider->setPageStep(10);
    // add them to layout
    QHBoxLayout *viewHLayout = new QHBoxLayout;
    QVBoxLayout *leftSideCanvas = new QVBoxLayout;
    leftSideCanvas->addWidget(glWidgetArea);
    leftSideCanvas->addWidget(timeSlider);
    leftSideCanvas->setContentsMargins(0,0,0,0);
    viewHLayout->addLayout(leftSideCanvas);
    /*************** control panel on the right *******************/
    contrastScrollBar = new QScrollBar(Qt::Orientation::Vertical);
    contrastScrollBar->setFocusPolicy(Qt::StrongFocus);
    contrastScrollBar->setRange(-100, 100);
    contrastScrollBar->setSingleStep(1);
    contrastScrollBar->setPageStep(10);
    QLabel *contrastSlider_Label = new QLabel("Contrast");
//    contrastSlider_Label->setWordWrap(true);
//    contrastSlider_Label->setAlignment(Qt::AlignTop);
//    QString s = "Contrast";
//    contrastSlider_Label->setText(s.split("", QString::SkipEmptyParts).join("\n"));
    QVBoxLayout *rightSideControlLayout = new QVBoxLayout;
    //rightSideControlLayout->addWidget(contrastSlider_Label);
    rightSideControlLayout->addWidget(contrastScrollBar, Qt::AlignCenter);

    //rightSideControlLayout->setContentsMargins(0,0,0,0);
    viewHLayout->addLayout(rightSideControlLayout);
    // Put the layout to the mainwindow
    grpBox4display_canvas->setLayout(viewHLayout);
}
void MainWindow::updateControlPanel(){
    timeSlider->setMinimum(0);
    long st = data4test->image4d->getTDim();
    if (st < 1) st = 1;
    timeSlider->setMaximum(st - 1);// start from 0
}
/**
 * brief: connect events
 * */
void MainWindow::connectSignal()
{
    if (widget_type != raycast_type)	return;

    if(importImageSeriesAct){
        connect(importImageSeriesAct, SIGNAL(triggered()), this, SLOT(importImageSeries()));
    }
    if (resetViewPoint){
        connect(resetViewPoint, SIGNAL(triggered()), glWidget_raycast, SLOT(setLightPositionZero()));
    }
    if (timeSlider) {
        connect(glWidget_raycast, SIGNAL(changeVolumeTimePoint(int)), timeSlider, SLOT(setValue(int)));
        connect(timeSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setVolumeTimePoint(int)));
    }
    if (contrastScrollBar) {
        connect(contrastScrollBar, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setContrast(int)));
    }
    connect(this, SIGNAL(signalDataLoaded()), this, SLOT(updateControlPanel())); // simply for easy reading
}
void MainWindow::importImageSeries()
{
    //QString filename = QString("/home/ccw/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/crop_embryo_data_500x500x30x40/8bits/embryo_TM481.tif");
    QString filename = QFileDialog::getOpenFileName(this);
    if (!filename.isEmpty()) {
        try
        {
            if(data4test->importData(filename)) //load data from given paths
            {
                emit signalDataLoaded();
            }
            // display in glWidget
            if (widget_type == my_simple_test_type){
            }
            else if (widget_type == raycast_type){
                glWidget_raycast->initVolume(data4test);
            }
            else //vaa3d_type
            {
                //glWidget = new EmT_GLWidget(data4test, this, filename);
            }

        }
        catch (...)
        {
            v3d_msg("You fail to import the specified image(s). The file may have certain problem, or is simply too big but you don't have enough memory.");
        }
    }
}
//void MainWindow::zeroViewPiont(){
//    if (widget_type == my_simple_test_type){
////        glWidget_simple->xRotationChanged(0);
////        glWidget_simple->yRotationChanged(0);
////        glWidget_simple->zRotationChanged(0);
//    }
//    else if (widget_type == raycast_type){
//        glWidget_raycast->set
//        //glWidget_raycast->xRotationChanged(0);
//        //glWidget_raycast->yRotationChanged(0);
//        //glWidget_raycast->zRotationChanged(0);
//    }
//    else //vaa3d_type
//    {
//        //
//    }
//}
MainWindow::~MainWindow()
{
}

