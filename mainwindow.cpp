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
    bndAxesShow = new QAction(tr("&BoundingBox"), this);
    //resetViewPoint->setShortcut(tr("Ctrl+R"));
    // status has not been defined
    //importImageFileAct->setStatusTip(tr("Import general image series"));
    editMenu->addAction(resetViewPoint);
    editMenu->addAction(bndAxesShow);
    /** ************* process menu ********************/
    processMenu = new QMenu(tr("&Process"), this);
    menuBar->addMenu(processMenu);
    segmentCell3d = new QAction(tr("&Cell Segmentation"), this);
    trackCell3d = new QAction(tr("&Cell Tracking"), this);
    processMenu->addAction(segmentCell3d);
    processMenu->addAction(trackCell3d);
    /** *************** Debug algorithms *****************/
    debugButton = new QAction(tr("Debug"), this);
    menuBar->addAction(debugButton);
    /** *************** About *************************/
    QMenu * aboutMenu = new QMenu(tr("About"), this);
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
    if (bndAxesShow){
        connect(bndAxesShow, SIGNAL(triggered()), glWidget_raycast, SLOT(setBnfAxesOnOff()));
    }
    if (timeSlider) {
        connect(glWidget_raycast, SIGNAL(changeVolumeTimePoint(int)), timeSlider, SLOT(setValue(int)));
        connect(timeSlider, SIGNAL(valueChanged(int)), this, SLOT(transferRGBAVolume(int)));
        //connect(timeSlider, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setVolumeTimePoint(int)));
    }
    if (contrastScrollBar) {
        connect(contrastScrollBar, SIGNAL(valueChanged(int)), glWidget_raycast, SLOT(setContrast(int)));
    }
    connect(this, SIGNAL(signalDataLoaded()), this, SLOT(updateControlPanel())); // simply for easy reading

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
void MainWindow::debugAlgorithm()
{
    if(cellTracker != nullptr){ // already finish the tracking
        tracking_result_exist = !tracking_result_exist;
        return;
    }
    //float kk = chi2inv(0.5, 2);
    algorithmDebug = true;
    this->data4test->debugMode = true;
    this->importImageSeries();
    //// send data to do segmentation on all frames
    for(int i = 0; i < glWidget_raycast->bufSize[4]; i++){
        glWidget_raycast->setVolumeTimePoint(i);
        this->sendData4Segment();
    }
    glWidget_raycast->setVolumeTimePoint(0);
    //// send segmentation results for cell linking
    this->sendData4Track();
}
void MainWindow::sendData4Segment()
{
    if(!data4test || !glWidget_raycast){
        QMessageBox::critical(0, "ASSERT", tr("data has not been imported or displayed"));
        return;
    }
    if(!cellSegmenter){
        cellSegmenter = new cellSegmentMain((void *)data4test->image4d->getRawData(),
                                            data4test->image4d->getDatatype(),
                                            glWidget_raycast->bufSize);
    }
    /// way 1: directly detect cells on the original data
    //cellSegmenter->processSingleFrameAndReturn(glWidget_raycast);
    /// way 2: try to load saved data. Detect cells if failed.
    cellSegmenter->processSingleFrameAndReturn(glWidget_raycast,
                   data4test->filelist.at(glWidget_raycast->curr_timePoint_in_canvas));
}

void MainWindow::sendData4Track()
{
    if(cellTracker == nullptr){
        chrono::steady_clock::time_point begin = chrono::steady_clock::now();

        cellTracker = new cellTrackingMain(*cellSegmenter);

        chrono::steady_clock::time_point end = chrono::steady_clock::now();
        qInfo("----------------time used: %.3f s", ((float)chrono::duration_cast<chrono::milliseconds>(end - begin).count())/1000);

        // label that the tracking results is OK for illustration
        tracking_result_exist = cellTracker->tracking_sucess;
        //transferRGBAVolume(0);
    }
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
