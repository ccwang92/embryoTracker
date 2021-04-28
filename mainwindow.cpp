#include "mainwindow.h"
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_message.h"
#include <chrono> // time elapsed
/**
 * brief: create main window
 * */
MainWindow::MainWindow()
{
    //QVBoxLayout *mainLayout = new QVBoxLayout();
    createControlWidgets();
    // initialize data4test as null
    data4test = new DataImporter(algorithmDebug); //do we need?
    setMenuBar(menuBar);
    setCentralWidget(grpBox4display_canvas);
    //layout->addWidget(rightSideControlLayout);
    //setLayout(mainLayout);
    this->setWindowTitle("BIQAV");
    connectSignal();
    cellSegmenter = 0;
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
/** ********** file menu define *****************/
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
/** ********** edit menu define *****************/
    editMenu = new QMenu(tr("&Edit"), this);
    menuBar->addMenu(editMenu);
    // reset view point
    resetViewPoint = new QAction(tr("&Reset ViewPoint"), this);
    resetViewPoint->setShortcut(tr("Ctrl+R"));
    //bndAxesShow = new QAction(tr("&Axes"), this);
    editMenu->addAction(resetViewPoint);
    //editMenu->addAction(bndAxesShow);
/** ************* process menu ********************/
    processMenu = new QMenu(tr("&Process"), this);
    menuBar->addMenu(processMenu);
    segmentCell3d = new QAction(tr("&Cell Segmentation"), this);
    trackCell3d = new QAction(tr("&Cell Tracking"), this);
    processMenu->addAction(segmentCell3d);
    processMenu->addAction(trackCell3d);
/** *************** Debug algorithms *****************/
    debugButton = new QAction(tr("Debug"), this);
    //menuBar->addAction(debugButton);
    processMenu->addAction(debugButton);
    /** *************** About *************************/
    QMenu * aboutMenu = new QMenu(tr("About"), this);
    menuBar->addMenu(aboutMenu);
/** ************* display grid *******************/
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
/** ************* control panel on the right *******************/
    QVBoxLayout *rightSideControlLayout = new QVBoxLayout;

    QWidget *volDisplayOptGroup = new QGroupBox("Control");
    //volDisplayOptGroup
    QGridLayout *layout_mainDisplayOptGroup = new QGridLayout(volDisplayOptGroup);
    contrastScrollBar = createContrastSlider();
    thresholdScrollBar = createThreshodSlider();
    axesCheckBox = new QCheckBox("Axes", volDisplayOptGroup);
    //bndboxCheckBox = new QCheckBox("Bounding box", volDisplayOptGroup);
    backgroundColorButton = new QPushButton("Change Canvas Color", volDisplayOptGroup);
    layout_mainDisplayOptGroup->addWidget(new QLabel("Threshold"), 1, 0, 1, 4);
    layout_mainDisplayOptGroup->addWidget(thresholdScrollBar, 1, 4, 1, 21-4);
    layout_mainDisplayOptGroup->addWidget(new QLabel("Contrast"), 2, 0, 1, 4);
    layout_mainDisplayOptGroup->addWidget(contrastScrollBar, 2, 4, 1, 21-4);
    layout_mainDisplayOptGroup->addWidget(axesCheckBox, 3, 0, 1, 3);
    //layout_mainDisplayOptGroup->addWidget(bndboxCheckBox, 3, 4, 1, 17);
    layout_mainDisplayOptGroup->addWidget(backgroundColorButton, 3, 3, 1, 18);
    rightSideControlLayout->addWidget(volDisplayOptGroup, Qt::AlignTop);
/** ************* Volume control panel on the right *******************/
    QWidget *volCutGroup = new QGroupBox("Volume Cut");
    QGridLayout *vol_yxzCntGroup = new QGridLayout(volCutGroup);
    int d1, d2, d3;
    d1 = 100;
    d2 = 100;
    d3 = 100;
    xcminSlider = createCutPlaneSlider(d1);
    xcmaxSlider = createCutPlaneSlider(d1); xcmaxSlider->setValue(d1);
    ycminSlider = createCutPlaneSlider(d2);
    ycmaxSlider = createCutPlaneSlider(d2); ycmaxSlider->setValue(d1);
    zcminSlider = createCutPlaneSlider(d3);
    zcmaxSlider = createCutPlaneSlider(d3); zcmaxSlider->setValue(d1);

    vol_yxzCntGroup->addWidget(new QLabel("X-cut"), 1, 0, 2, 4);
    vol_yxzCntGroup->addWidget(xcminSlider, 1, 4, 1, 17);
    vol_yxzCntGroup->addWidget(xcmaxSlider, 2, 4, 1, 17);

    vol_yxzCntGroup->addWidget(new QLabel("Y-cut"), 3, 0, 2, 4);
    vol_yxzCntGroup->addWidget(ycminSlider, 3, 4, 1, 17);
    vol_yxzCntGroup->addWidget(ycmaxSlider, 4, 4, 1, 17);

    vol_yxzCntGroup->addWidget(new QLabel("Z-cut"), 5, 0, 2, 4);
    vol_yxzCntGroup->addWidget(zcminSlider, 5, 4, 1, 17);
    vol_yxzCntGroup->addWidget(zcmaxSlider, 6, 4, 1, 17);

    rightSideControlLayout->addWidget(volCutGroup, Qt::AlignTop);
/** ************* Save Path panel on the right *******************/
    QWidget *saveResultGroup = new QGroupBox("Save");
    QGridLayout *saveResultsOptGroup = new QGridLayout(saveResultGroup);
    saveSegCheckBox = new QCheckBox("Save Segmentation Results", saveResultGroup);
    saveTrackCheckBox = new QCheckBox("Save Tracking Results", saveResultGroup);
    changeSaveFolderButton = new QPushButton("Change Folder", saveResultGroup);
    saveFolder = new QLineEdit;
    saveFolderName = "";
    saveFolder->setText(saveFolderName);
    saveResultsOptGroup->addWidget(saveSegCheckBox, 1, 0, 1, 21);
    saveResultsOptGroup->addWidget(saveTrackCheckBox, 2, 0, 1, 21);
    saveResultsOptGroup->addWidget(new QLabel("Save Folder"), 3, 0, 1, 8);
    saveResultsOptGroup->addWidget(saveFolder, 3, 8, 1, 13);
    saveResultsOptGroup->addWidget(changeSaveFolderButton, 4, 10, 1, 11);
    rightSideControlLayout->addWidget(saveResultGroup, Qt::AlignTop);
// --------------------------------------------------------------------//
    viewHLayout->addLayout(rightSideControlLayout);
    // Put the layout to the mainwindow
    grpBox4display_canvas->setLayout(viewHLayout);
}

QAbstractSlider *MainWindow::createCutPlaneSlider(int maxval, Qt::Orientation hv)
{
//    QSlider *slider = new QSlider(hv);
//    slider->setRange(0, maxval);
//    slider->setSingleStep(1);
//    slider->setPageStep(10);
//    slider->setTickInterval(10);
//    slider->setTickPosition(QSlider::TicksRight);
    QScrollBar *slider = new QScrollBar(hv);
    slider->setRange(0, maxval);
    slider->setSingleStep(1);
    slider->setPageStep(10);

    //slider->setValue(0);
    return slider;
}
QScrollBar* MainWindow::createContrastSlider(){
    QScrollBar *cstScrollBar = new QScrollBar(Qt::Orientation::Horizontal);
    //cstScrollBar->setFocusPolicy(Qt::StrongFocus);
    cstScrollBar->setRange(-100, 100);
    cstScrollBar->setSingleStep(1);
    cstScrollBar->setPageStep(10);
    return cstScrollBar;
}

QScrollBar* MainWindow::createThreshodSlider(){
    QScrollBar *thresScrollBar = new QScrollBar(Qt::Orientation::Horizontal);
    thresScrollBar->setRange(0, 100);
    thresScrollBar->setSingleStep(1);
    thresScrollBar->setPageStep(10);
    return thresScrollBar;
}
//QCheckBox* createAxesCheckBox(){
//    QCheckBox *axesCkBox = new QCheckBox();
//    return axesCkBox;
//}
//QCheckBox* createBndBoxCheckBox(){
//    QCheckBox *bndBoxCkBox = new QCheckBox();
//    return bndBoxCkBox;
//}
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
//    if (bndAxesShow){
//        connect(bndAxesShow, SIGNAL(triggered()), glWidget_raycast, SLOT(setBnfAxesOnOff()));
//        //connect(bndAxesShow, SIGNAL(triggered()), axesCheckBox, SLOT(toggle()));
//    }
    if (timeSlider) {
        connect(glWidget_raycast, SIGNAL(changeVolumeTimePoint(int)), timeSlider, SLOT(setValue(int)));
        //connect(timeSlider, SIGNAL(valueChanged(int)), this, SLOT(transferRGBAVolume(int)));
        connect(timeSlider, SIGNAL(valueChanged(int)), this, SLOT(cellTraceAppend(int)));
        //connect(timeSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setVolumeTimePoint(int)));
    }
    if (contrastScrollBar) {
        connect(contrastScrollBar, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setContrast(int)));
    }
    if (thresholdScrollBar) {
        connect(thresholdScrollBar, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setThreshold(int)));
    }
//    xcminSlider = createCutPlaneSlider(d1);
//    xcmaxSlider = createCutPlaneSlider(d1);
//    ycminSlider = createCutPlaneSlider(d2);
//    ycmaxSlider = createCutPlaneSlider(d2);
//    zcminSlider = createCutPlaneSlider(d3);
//    zcmaxSlider = createCutPlaneSlider(d3);
    connect(xcminSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeXMIN(int)));
    connect(xcmaxSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeXMAX(int)));
    connect(ycminSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeYMIN(int)));
    connect(ycmaxSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeYMAX(int)));
    connect(zcminSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeZMIN(int)));
    connect(zcmaxSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setRangeZMAX(int)));

    connect(this, SIGNAL(signalDataLoaded()), this, SLOT(updateControlPanel())); // simply for easy reading

    if(backgroundColorButton){
        connect(backgroundColorButton, SIGNAL(clicked()), this, SLOT(setBackgroundColor()));
    }
    if (axesCheckBox) {
        connect(axesCheckBox, SIGNAL(toggled(bool)), glWidget_raycast, SLOT(setBnfAxesOnOff()));
        //connect(glWidget_raycast, SIGNAL(changeBnfAxesOnOff(bool)), this, SLOT(setAxesCheckBox(bool)));
    }
    /** cell segmentation and tracking algorithm call ***/
    if (segmentCell3d){
        connect(segmentCell3d, SIGNAL(triggered()), this, SLOT(sendData4Segment()));
    }
    if (trackCell3d){
        connect(trackCell3d, SIGNAL(triggered()), this, SLOT(sendData4Track()));
    }
    if (debugButton){
        connect(debugButton, SIGNAL(triggered()), this, SLOT(debugAlgorithm()));
    }
}
void MainWindow::importImageSeries()
{
    QString filename;
    if (algorithmDebug){
        filename = debugDataPath;
    }else{
        filename = QFileDialog::getOpenFileName(this);
    }
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
void MainWindow::debugBatchFusion(){
    // directly work on 4 folders
    QString dataFolderName =  QFileDialog::getExistingDirectory(this, tr("Open Image Directory"), "/home",
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    QString resFolderName =  QFileDialog::getExistingDirectory(this, tr("Open Batch Results Directory"), dataFolderName,
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    
    qFatal("Debug successed! No need to continue");
}
void MainWindow::debugAlgorithm()
{
    // directly work on 4 folders
    QString filename = QString::fromStdString("/home/ccw/storage/debug_batch_merge/TM450_489/backleft/1.tif");
    //QString filename = QFileDialog::getOpenFileName(this);
    int x = filename.lastIndexOf('/');
    QString folderName = filename.left(x);
    x = folderName.lastIndexOf('/');
    folderName = folderName.left(x);
    QDirIterator it(folderName, QDir::Dirs | QDir::NoDotAndDotDot); //QDirIterator::Subdirectories);//QStringList() << "*.jpg", QDir::Files,
    while (it.hasNext()){
        debugDataPath = it.next();
        qInfo()<<debugDataPath;
        if(!debugDataPath.isEmpty()){
            QFileInfo check_file(QDir(debugDataPath).filePath("movieInfo.txt"));
            if(check_file.exists() && check_file.isFile()){
                continue;
            }
            debugDataPath = QDir(debugDataPath).filePath("1.tif");
            algorithmDebug = true;
            this->data4test->debugMode = true;
            //QString fileName =
            this->importImageSeries();
            //// send data to do segmentation on all frames
            for(int i = 0; i < glWidget_raycast->bufSize[4]; i++){
                glWidget_raycast->curr_timePoint_in_canvas = i;
                //// segement
                this->sendData4Segment();
                qInfo("The #%d/%ld frame are finished!", i, glWidget_raycast->bufSize[4]);
            }
        //    glWidget_raycast->setVolumeTimePoint(0);
            //// send segmentation results for cell linking
        //    this->sendData4Track();
            qInfo("start tracking");
            chrono::steady_clock::time_point begin = chrono::steady_clock::now();
            cellTracker = new cellTrackingMain(*cellSegmenter, data4test->filelist);
            chrono::steady_clock::time_point end = chrono::steady_clock::now();
            qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);
            delete cellSegmenter;
            cellSegmenter = nullptr;
            delete cellTracker;
            cellTracker = nullptr;
        }else{

        }
    }
    qFatal("Debug successed! No need to continue");
}
void MainWindow::sendData4Segment()
{
    if(!data4test->image4d || !glWidget_raycast->getDataImporter()){
        QMessageBox::critical(0, "ASSERT", tr("data has not been imported or displayed"));
        return;
    }
    if(cellSegmenter == nullptr){
        cellSegmenter = new cellSegmentMain((void *)data4test->image4d->getRawData(),
                                            data4test->image4d->getDatatype(),
                                            glWidget_raycast->bufSize);
    }
    /// way 1: directly detect cells on the original data
    //cellSegmenter->processSingleFrameAndReturn(glWidget_raycast);
    /// way 2: try to load saved data. Detect cells if failed.
    cellSegmenter->processSingleFrameAndReturn(glWidget_raycast,
                   data4test->filelist.at(glWidget_raycast->curr_timePoint_in_canvas), seg4track);
    //// display results in canvas
    int i = glWidget_raycast->curr_timePoint_in_canvas;
    long sz_single_frame = cellSegmenter->data_rows_cols_slices[0]*cellSegmenter->data_rows_cols_slices[1]*
            cellSegmenter->data_rows_cols_slices[2];
    unsigned char *ind = (unsigned char*)cellSegmenter->normalized_data4d.data + sz_single_frame*i; // sub-matrix pointer
    Mat *single_frame = new Mat(3, cellSegmenter->normalized_data4d.size, CV_8U, ind);
    Mat4b rgb_mat4display;
    label2rgb3d(cellSegmenter->cell_label_maps[i], *single_frame, rgb_mat4display);
    glWidget_raycast->setMode("Alpha blending rgba");
    glWidget_raycast->getRenderer()->transfer_volume((unsigned char *)rgb_mat4display.data, 0, 255, cellSegmenter->data_rows_cols_slices[1],
            cellSegmenter->data_rows_cols_slices[0], cellSegmenter->data_rows_cols_slices[2], 4);

}

void MainWindow::sendData4Track()
{
    if(!data4test->image4d || !glWidget_raycast->getDataImporter()){
        QMessageBox::critical(0, "ASSERT", tr("data has not been imported or displayed"));
        return;
    }
    if(cellTracker == nullptr){
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();
        QString fileName = data4test->filelist.at(0);
        QString movieInfo_txt_name = fileName.left(fileName.lastIndexOf('/')) + "/movieInfo.txt";
        QFileInfo check_file(movieInfo_txt_name);
        if(check_file.exists() && check_file.isFile()){
            if(!cellSegmenter){ // already finish the tracking
                cellSegmenter = new cellSegmentMain((void *)data4test->image4d->getRawData(),
                                                    data4test->image4d->getDatatype(),
                                                    glWidget_raycast->bufSize);
            }
            vector<int> yxzt_sz;
            yxzt_sz.reserve(4);
            yxzt_sz.push_back(glWidget_raycast->bufSize[1]);
            yxzt_sz.push_back(glWidget_raycast->bufSize[0]);
            yxzt_sz.push_back(glWidget_raycast->bufSize[2]);
            yxzt_sz.push_back(glWidget_raycast->bufSize[4]);
            //glWidget_raycast->setVolumeTimePoint(0);
            cellTracker = new cellTrackingMain(yxzt_sz, data4test->filelist);
        }else{
            if(!cellSegmenter){ // already finish the tracking
                //// send data to do segmentation on all frames
                seg4track = true;
                for(int i = 0; i < glWidget_raycast->bufSize[4]; i++){
                    //glWidget_raycast->setVolumeTimePoint(i);
                    glWidget_raycast->curr_timePoint_in_canvas = i;
                    this->sendData4Segment();
                    qInfo("The #%d/%ld frame are finished!", i, glWidget_raycast->bufSize[4]);
                }
                seg4track = false;
            }
            //glWidget_raycast->setVolumeTimePoint(0);
            cellTracker = new cellTrackingMain(*cellSegmenter, data4test->filelist);
        }
        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);

        // label that the tracking results is OK for illustration
        tracking_result_exist = cellTracker->tracking_sucess;
    }else{
        tracking_result_exist = !tracking_result_exist;
    }
    //transferRGBAVolume(0);
    cellTraceAppend(0);
}
void MainWindow::cellTraceAppend(int t){
    glWidget_raycast->show_track_result = tracking_result_exist;
    if(tracking_result_exist){ // transfer the volume to glWidget_raycast->rgb_frame
        vector<int> yxzt_sz;
        yxzt_sz.reserve(4);
        yxzt_sz.push_back(glWidget_raycast->bufSize[1]);
        yxzt_sz.push_back(glWidget_raycast->bufSize[0]);
        yxzt_sz.push_back(glWidget_raycast->bufSize[2]);
        yxzt_sz.push_back(glWidget_raycast->bufSize[4]);
        glWidget_raycast->rgb_frame = Mat(); // clear the content by assign an empty mat
        // re-set the rgb_frame
        long sz_single_frame = cellSegmenter->data_rows_cols_slices[0]*
                cellSegmenter->data_rows_cols_slices[1]*cellSegmenter->data_rows_cols_slices[2];
        unsigned char *ind = (unsigned char*)cellSegmenter->normalized_data4d.data + sz_single_frame*t; // sub-matrix pointer
        Mat *single_frame = new Mat(3, cellSegmenter->normalized_data4d.size, CV_8U, ind);
        //label2rgb3d(cellSegmenter->cell_label_maps[t], *single_frame, glWidget_raycast->rgb_frame);

        // build a map based on tracking results
        // 1. build color map and extract trace locations
        if(colormap4tracking_res.empty()){
            //int max_cell_num = vec_max(cellSegmenter->number_cells);
            int color_num =  cellTracker->movieInfo.tracks.size(); //MAX((size_t)max_cell_num, cellTracker->movieInfo.tracks.size());
            colorMapGen((double)color_num, colormap4tracking_res);
            cellTracker->extractTraceLocations(yxzt_sz);
        }
//        label2rgb3d(cellSegmenter->cell_label_maps[t], *single_frame, colormap4tracking_res, glWidget_raycast->rgb_frame);
        // 2. build label map
        Mat1i mapedLabelMap = Mat::zeros(3, yxzt_sz.data(), CV_32S);

        for(int j=0; j<cellTracker->movieInfo.tracks.size(); j++){
            if(cellTracker->movieInfo.tracks[j].size()<=5){ // if the trace has stopped before time t
                continue;
            }
            int end_time = cellTracker->movieInfo.frames[*cellTracker->movieInfo.tracks[j].rbegin()];
            int start_time = cellTracker->movieInfo.frames[*cellTracker->movieInfo.tracks[j].begin()];
            if(end_time < t || start_time > t){ // if the trace has stopped before time t
                continue;
            }
            for(int i=0; i<=t; i++){
                for(auto idx : cellTracker->trace_sets[i][j]){
                    mapedLabelMap.at<int>(idx) = j + 1;
                }
            }
        }
        label2rgb3d(mapedLabelMap, *single_frame, colormap4tracking_res, glWidget_raycast->rgb_frame);
    }
    //qInfo("3-Done!");
    glWidget_raycast->setVolumeTimePoint(t);
}
void MainWindow::transferRGBAVolume(int t){
    glWidget_raycast->show_track_result = tracking_result_exist;
    if(tracking_result_exist){ // transfer the volume to glWidget_raycast->rgb_frame
        glWidget_raycast->rgb_frame = Mat(); // clear the content by assign an empty mat
        // re-set the rgb_frame
        long sz_single_frame = cellSegmenter->data_rows_cols_slices[0]*
                cellSegmenter->data_rows_cols_slices[1]*cellSegmenter->data_rows_cols_slices[2];
        unsigned char *ind = (unsigned char*)cellSegmenter->normalized_data4d.data + sz_single_frame*t; // sub-matrix pointer
        Mat *single_frame = new Mat(3, cellSegmenter->normalized_data4d.size, CV_8U, ind);
        //label2rgb3d(cellSegmenter->cell_label_maps[t], *single_frame, glWidget_raycast->rgb_frame);

        // build a map based on tracking results
        // 1. build color map
        if(colormap4tracking_res.empty()){
            int max_cell_num = vec_max(cellSegmenter->number_cells);
            int color_num = MAX((size_t)max_cell_num, cellTracker->movieInfo.tracks.size());
            colorMapGen((double)color_num, colormap4tracking_res);
        }
//        label2rgb3d(cellSegmenter->cell_label_maps[t], *single_frame, colormap4tracking_res, glWidget_raycast->rgb_frame);
        // 2. build label map
        Mat1i mapedLabelMap = Mat::zeros(cellSegmenter->cell_label_maps[t].dims,
                                         cellSegmenter->cell_label_maps[t].size, CV_32S);
        vector<bool> color_used (colormap4tracking_res.size[0] + 1, false);
        //qDebug("%d-%d", colormap4tracking_res.size[0], colormap4tracking_res.size[1]);
        FOREACH_i_MAT(cellSegmenter->cell_label_maps[t]){
            size_t idx = cellSegmenter->cell_label_maps[t].at<int>(i);
            if(idx == 0) continue;
            idx --;
            idx += t == 0 ? 0 : cellTracker->cumulative_cell_nums[t-1];
            if(cellTracker->movieInfo.nodes[idx].nodeId2trackId >= 0){
                color_used[cellTracker->movieInfo.nodes[idx].nodeId2trackId] = true;
            }
        }
        size_t next_available = 0;
        vector<int> mapped_idx (color_used.size(), -1);
        FOREACH_i_MAT(cellSegmenter->cell_label_maps[t]){
            size_t org_idx = cellSegmenter->cell_label_maps[t].at<int>(i);
            if(org_idx == 0) continue;
            size_t idx = org_idx-1;
            idx += t == 0 ? 0 : cellTracker->cumulative_cell_nums[t-1];
            int track_id = cellTracker->movieInfo.nodes[idx].nodeId2trackId;
            if(track_id >= 0 && cellTracker->movieInfo.tracks[track_id].size() > 1){
                mapedLabelMap.at<int>(i) = track_id + 1;
            }else{
                continue;
                //                if(mapped_idx[org_idx] >= 0){
                //                    mapedLabelMap.at<int>(i) = mapped_idx[org_idx];
                //                }else{
                //                    while(color_used[next_available]){
                //                        next_available ++;
                //                    }
                //                    mapedLabelMap.at<int>(i) = next_available + 1;
                //                    mapped_idx[org_idx] = next_available + 1;
                //                    if(next_available >= color_used.size()){
                //                        qDebug("%d", next_available);
                //                    }
                //                    color_used[next_available] = true;
                //                }
            }
        }
        label2rgb3d(mapedLabelMap, *single_frame, colormap4tracking_res, glWidget_raycast->rgb_frame);
//        label2rgb3d(mapedLabelMap, *single_frame, glWidget_raycast->rgb_frame);
//        glWidget_raycast->setMode("Alpha blending rgba");
//        glWidget->getRenderer()->transfer_volume((unsigned char *)rgb_mat4display.data, 0, 255, data_rows_cols_slices[1],
//                data_rows_cols_slices[0], data_rows_cols_slices[2], 4);
    }
    glWidget_raycast->setVolumeTimePoint(t);
}

void MainWindow::setBackgroundColor(){
    const QColor colour = QColorDialog::getColor(glWidget_raycast->getBackground(), this, "Select background colour");

    if (colour.isValid()) {
        glWidget_raycast->setBackground(colour);
    }
}

//void MainWindow::setAxesCheckBox(bool axesOn){
//    axesCheckBox->setChecked(axesOn);
//}
