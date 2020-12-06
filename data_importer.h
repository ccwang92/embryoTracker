#ifndef DataImporter_H
#define DataImporter_H
#include <QtGui>
//#include "my4dimage.h"
#include "src_3rd/basic_c_fun/basic_4dimage.h"
//#include "src_3rd/basic_c_fun/basic_landmark.h"


class DataImporter
{
public:
    DataImporter(){};
    ~DataImporter(){};

    void importData(QString filename);
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
public:
    QString openFileNameLabel = QString("");;
    Image4DSimple* image4d = 0;
};

#endif // DataImporter_H
