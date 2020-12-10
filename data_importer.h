#ifndef DataImporter_H
#define DataImporter_H
#include <QtGui>
//#include "my4dimage.h"
#include "src_3rd/basic_c_fun/basic_4dimage.h"
//#include "src_3rd/basic_c_fun/basic_landmark.h"


class DataImporter//:public QObject //for emmit signal,but seems no need
{
public:
    DataImporter(){};
    ~DataImporter(){};

    bool importData(QString filename);
    bool importGeneralImgSeries(
            const QStringList & mylist, TimePackType timepacktype);

    bool importGeneralImageFile(QString filename);
    QStringList importSeriesFileList_addnumbersort(
            const QString & individualFileName,
            TimePackType & timepacktype);
    bool readSingleImageFile(
            char *imgSrcFile, unsigned char * & data1d,
            V3DLONG * & sz, ImagePixelType & datatype);
    void cleanData();
signals:
    //Nope

public:
    QString openFileNameLabel = QString("");
    Image4DSimple* image4d = 0;

    double * p_vmax = NULL, * p_vmin = NULL; //whole volume max/min values. Use pointer to handle multiple channels separately
    void updateMinMax(V3DLONG nFrame);
    bool updateminmaxvalues();
    //double getChannalMinIntensity(V3DLONG channo);
    //double getChannalMaxIntensity(V3DLONG channo);

};

#endif // DataImporter_H
