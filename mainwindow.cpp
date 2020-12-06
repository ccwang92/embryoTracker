#include "mainwindow.h"
#include <QMenuBar>
#include <QFileDialog>
#include <QStatusBar>
#include "data_importer.h"
#include "src_3rd/basic_c_fun/v3d_message.h"

//#include "src_3rd/v3d/xformwidget.h"

MainWindow::MainWindow()
{
    /************ Menu define *****************/
    fileMenu = menuBar()->addMenu(tr("&File"));
    // import image stacks
    importImageSeriesAct = new QAction(tr("&Import time series to an image stack..."), this);
    importImageSeriesAct->setShortcut(tr("Ctrl+I"));
    // status has not been defined
    //importImageFileAct->setStatusTip(tr("Import general image series"));
    connect(importImageSeriesAct, SIGNAL(triggered()), this, SLOT(importImageSeries()));
    fileMenu->addAction(importImageSeriesAct);
    // separator
    fileMenu->addSeparator();
    // exit the program
    fileMenu->addAction(exitAct);

    // initialize data4test as null
    data4test = new DataImporter();
}
void MainWindow::importImageSeries()
{
    QString filename = QString("/home/ccw/Insync/ccwang@vt.edu/Google Drive/Projects/embyo_analysis/data/crop_embryo_data_500x500x30x40/8bits/embryo_TM481.tif");
    //QString filename = QFileDialog::getOpenFileName(this);
    if (!filename.isEmpty()) {
        try
        {
            data4test->importData(filename);
        }
        catch (...)
        {
            v3d_msg("You fail to import the specified image(s). The file may have certain problem, or is simply too big but you don't have enough memory.");
        }
    }
}
MainWindow::~MainWindow()
{
}

