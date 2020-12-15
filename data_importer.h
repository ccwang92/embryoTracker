#ifndef DataImporter_H
#define DataImporter_H
#include <QtGui>
//#include "my4dimage.h"
#include "src_3rd/basic_c_fun/basic_4dimage.h"
//#include "src_3rd/basic_c_fun/basic_landmark.h"
#include "src_3rd/basic_c_fun/color_xyz.h"
#include <QMessageBox>
#define MESSAGE_ASSERT(s) \
{\
    if (!(s)) \
        QMessageBox::critical(0, "ASSERT", QObject::tr("ASSERT(%1) in file(%2) at line(%3)").arg(#s).arg(__FILE__).arg(__LINE__)); \
    Q_ASSERT(s); \
}

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
    // update the input gray image data into 8-bit color (RGBA) data
    void rgba3d_r2gray(RGBA8* rgbaBuf, V3DLONG bufSize[5]);
    void data4dp_to_rgba3d(Image4DProxy<Image4DSimple>& img4dp, V3DLONG dim5,
            V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
            V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
            RGBA8* rgbaBuf, V3DLONG bufSize[5]);
//    void data4dp_to_rgba3d(unsigned char* data4dp, V3DLONG dim1, V3DLONG dim2, V3DLONG dim3, V3DLONG dim4, V3DLONG dim5,
//            V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
//            V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
//            RGBA8* rgbaBuf, V3DLONG bufSize[5]);
    float sampling3dUINT8(Image4DProxy<Image4DSimple>& img4dp,
            V3DLONG c,
            V3DLONG x, V3DLONG y, V3DLONG z, V3DLONG dx, V3DLONG dy, V3DLONG dz);
    float sampling3dUINT8_2(Image4DProxy<Image4DSimple>& img4dp,
                    V3DLONG c,
                    V3DLONG x, V3DLONG y, V3DLONG z, V3DLONG dx, V3DLONG dy, V3DLONG dz, V3DLONG dxyz);

    //double getChannalMinIntensity(V3DLONG channo);
    //double getChannalMaxIntensity(V3DLONG channo);

};

#endif // DataImporter_H
