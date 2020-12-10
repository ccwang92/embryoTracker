#include "data_importer.h"
#include <stdio.h>
#include "src_3rd/v3d/import_images_tool_dialog.h"

#include "src_3rd/basic_c_fun/stackutil.h"
//#include "src_3rd/basic_c_fun/volimg_proc_declare.h"
#include "src_3rd/io/io_bioformats.h"
#include "src_3rd/basic_c_fun/v3d_message.h"
//#include "v3d_message.h"


bool DataImporter::importData(QString filename)
{
    if (!filename.isEmpty())
    {
        openFileNameLabel = filename;
        TimePackType timepacktype;
        QStringList mylist = importSeriesFileList_addnumbersort(filename, timepacktype);

        if (importGeneralImgSeries(mylist, timepacktype))
        {
            //updateMinMax(0); // use the first frame to update the minmax intensity value
            v3d_msg("File loaded!");
            return true;
        }
        else
            v3d_msg("Fail loading!");
    }
    else
        v3d_msg("Error file name!");
    return false;
}

bool DataImporter::importGeneralImgSeries(const QStringList & mylist, TimePackType timepacktype)
{
    //foreach (QString qs, myList)  qDebug() << qs;
    V3DLONG ntime = mylist.size();
    if (ntime <1)
    {
        v3d_msg("The import list is empty. do nothing.\n");
        return false;
    }

    //if there are files to import, then clean data

    if (image4d)
        cleanData();

    image4d = new Image4DSimple;
    if (!image4d)  return false;

    //now we can simply read files one by one, and arrange them in term of (1) the color channel and (2) z-planes. At this stage we need to verify each plane has the same size, and NOT colored image!!

    V3DLONG nsz0=0, nsz1=0;
    ImagePixelType cur_datatype=V3D_UNKNOWN, ndatatype=V3D_UNKNOWN;
    V3DLONG ncolors=0, nthick=0;
    V3DLONG pack_color=0, pack_z=0;

    for (V3DLONG i = 0; i < ntime; ++i)
    {
        QString tmpfileInfo = mylist.at(i);
        printf("importing %i file: {%s}\n", i, qPrintable(tmpfileInfo));

        unsigned char * cur_data1d=0;
        V3DLONG * cur_sz=0;
        if (!readSingleImageFile(tmpfileInfo.toUtf8().data(), cur_data1d, cur_sz, cur_datatype))
        {
            v3d_msg("Error occurs in reading the file. Exit importing.\n");
            return false;
        }
        if (!cur_data1d || !cur_sz)
        {
            v3d_msg("Error occurs in reading the file. Exit importing.\n");
            return false;
        }

        //-----------------------------------------------------------------------
        // 090731 RZC: (3D time series --> 4D color image) packed time in Z dim.
        //-----------------------------------------------------------------------

        if (i==0) ncolors = cur_sz[3];
        if (cur_sz[3]<=0 || cur_sz[3]!=ncolors)
        {
            printf("The current file has invalid or different colors [=%ld] from first section [=%ld]. Exit importing.\n", cur_sz[3], ncolors);
            v3d_msg("The current file has invalid or different colors\n");
            if (cur_data1d) {delete []cur_data1d; cur_data1d=0;}
            if (cur_sz) {delete []cur_sz; cur_sz=0;}
            return false;
        }
        if (i==0) ndatatype = cur_datatype;
        if (cur_datatype != ndatatype)
        {
            printf("The current file has different data type [=%ld] from first section [=%ld]. Exit importing.\n", cur_datatype, ndatatype);
            v3d_msg("The current file has different data type\n");
            if (cur_data1d) {delete []cur_data1d; cur_data1d=0;}
            if (cur_sz) {delete []cur_sz; cur_sz=0;}
            return false;
        }


        if (i==0)
        {
            nsz0 = cur_sz[0]; nsz1 = cur_sz[1]; nthick = cur_sz[2];
            if (timepacktype==TIME_PACK_Z)
            {
                pack_z     = nthick*ntime;
                pack_color = ncolors;
            }
            else // TIME_PACK_C
            {
                pack_z     = nthick;
                pack_color = ncolors*ntime;
            }

            if (image4d->createImage(nsz0, nsz1, pack_z, pack_color,
                                     cur_datatype)==false)
            {
                v3d_msg("Fail to allocate memory for the image stack. Exit importing.\n");
                return false;
            }
            image4d->setTDim( ntime );
            image4d->setTimePackType( timepacktype );
        }
        else
        {
            if (cur_sz[0]!=nsz0 || cur_sz[1]!=nsz1 || cur_sz[2]!=nthick)
            {
                printf("The current image has a different [width, height, thick]=[%ld, %ld, %ld] from the first section [%ld, %ld, %ld]. Exit importing.\n",
                       cur_sz[0], cur_sz[1], cur_sz[2], nsz0, nsz1, nthick);
                v3d_msg("The current image has a different [width, height, thick]\n");
                if (cur_data1d) {delete []cur_data1d; cur_data1d=0;}
                if (cur_sz) {delete []cur_sz; cur_sz=0;}
                return false;
            }
        }

        //now copy data of different planes into the 5D stack
        V3DLONG element_bytes = image4d->getUnitBytes();
        V3DLONG cur_time = i;
        V3DLONG block_size = (nthick*nsz0*nsz1)*(element_bytes);
        for (V3DLONG cur_ch=0; cur_ch<ncolors; cur_ch++)
        {
            unsigned char * cur_data1d_block = cur_data1d + (cur_ch)*block_size;
            unsigned char * cur_target1d_block;
            if (timepacktype==TIME_PACK_Z)
            {
                cur_target1d_block = image4d->getRawData() + (cur_ch*ntime + cur_time)*block_size;
            }
            else
            {
                cur_target1d_block = image4d->getRawData() + (cur_ch + cur_time*ncolors)*block_size;
            }

            //for (V3DLONG j=0; j<(block_size); j++)   cur_target1d_block[j] = cur_data1d_block[j];
            memcpy(cur_target1d_block, cur_data1d_block, block_size);
        }

        //now delete the temporary image data
        if (cur_data1d) {delete []cur_data1d; cur_data1d=0;}
        if (cur_sz) {delete []cur_sz; cur_sz=0;}

    }
    printf("Finished importing data. Now img data size = [%ld, %ld, %ld, %ld]\n", image4d->getXDim(), image4d->getYDim(), image4d->getZDim(), image4d->getCDim());
    if (image4d->getTDim()>1 && image4d->getTimePackType()==TIME_PACK_Z)
    {
        printf("Packed time point in [z = %ld * %ld]\n", image4d->getTDim(), image4d->getZDim()/image4d->getTDim());
    }
    if (image4d->getTDim()>1 && image4d->getTimePackType()==TIME_PACK_C)
    {
        printf("Packed time point in [c = %ld * %ld]\n", image4d->getTDim(), image4d->getCDim()/image4d->getTDim());
    }

    return true;
}


QStringList DataImporter::importSeriesFileList_addnumbersort(const QString & individualFileName, TimePackType & packtype)
{
    QStringList myList;
    myList.clear();

    //Get the image files namelist in the directory

    QFileInfo fileInfo(individualFileName);
    QString curFilePath = fileInfo.path();
        QString curSuffix = fileInfo.suffix();

    QDir dir(curFilePath);
    if (!dir.exists())
    {
        qWarning("Cannot find the directory");
        return myList;
    }

    QStringList imgfilters;
        imgfilters.append("*." + curSuffix);
    foreach (QString file, dir.entryList(imgfilters, QDir::Files, QDir::Name))
    {
        myList += QFileInfo(dir, file).absoluteFilePath();
    }

    //sort image sequence by numeric order instead of alphanumeric order
    //e.g. "a9.tiff" "a10.tiff" "b1.tiff"
    QStringList sortedList, tmpList;

    //-----------------------------------------------------------------------
    // 090731 RZC: fixed numerically sorting file names list, for DataImporter::importGeneralImgSeries
    //-----------------------------------------------------------------------
    QString fileNameStr, fileNameDigits;	//absolute file name is separated to 2 parts: strings and digits
        QRegExp r("(\\d+)");		//find digits
    QMap<V3DLONG, QString> mapList;

    mapList.clear();
    for(V3DLONG i=0; i<myList.size(); ++i)
    {
        fileNameStr = myList.at(i);
        QFileInfo fileFullName(myList.at(i));
        QString fileFullNameStr = fileFullName.completeBaseName();


        //extract the fileNameDigits from fileNameStr
        //e.g. "a9_b2009051801.tif.raw" into "a9_b.tif.raw" & "2009051801"

        V3DLONG pos = 0;
        fileNameDigits = "";
        while ((pos = r.indexIn(fileFullNameStr, pos)) != -1)
        {
                    fileNameDigits = r.cap(1);
                    pos += r.matchedLength();
        }

        if (fileNameDigits.isEmpty()) continue;


        V3DLONG num = fileNameDigits.toULong();
        mapList.insert(num, fileNameStr);
    }
    // must be sorted by QMap
    myList = mapList.values();
    //foreach (QString qs, myList)  qDebug() << qs;

    //no need to pop-out a dialog if no meaningful file has been detected. 131017
    if (myList.isEmpty())
    {
        v3d_msg("It seems no file contains a digit-portion in the file name. Naming convention should be sth like xxx_000001.yyy, xxx_000002.yyy, .... Check your data before importing.");
        return myList;
    }

    //Get the tiff image sequence by usr interaction

    ImportImgPara p;
    p.countImg = myList.size();

    import_images_tool_Dialog d(curFilePath);

    //need to update the information based on the current myList info. 131017
    d.numimgBox->setMaximum(p.countImg);
    d.numimgBox->setValue(p.countImg);
    d.numimgBox->setMinimum(p.countImg);

    d.startimgBox->setMaximum(p.countImg);
    d.startimgBox->setValue(1);

    d.incBox->setMaximum(p.countImg);

    d.endimgBox->setMaximum(p.countImg);
    d.endimgBox->setValue(p.countImg);

    int res = d.exec();
    if (res==QDialog::Accepted)
    {
        d.fetchData(&p);

        //get the QStringList
        myList = myList.filter(p.filt);

        tmpList.clear();
        for (V3DLONG i = p.startImg-1; i<= p.endImg-1; i+=p.inc)
            tmpList << myList.at(i);//.toLocal8Bit().constData();

        myList = tmpList;
        packtype = (p.packType==0)? TIME_PACK_Z : TIME_PACK_C;
    }
    else
    {
        myList.clear();
    }

    return myList;
}

bool DataImporter::readSingleImageFile(char *imgSrcFile, unsigned char * & data1d, V3DLONG * & sz, ImagePixelType & datatype)
{
    datatype = V3D_UNKNOWN;
    int dt = 0;
    if (loadImage(imgSrcFile, data1d, sz,  dt))
    {
        if (dt==1) datatype = V3D_UINT8;
        else if (dt==2) datatype = V3D_UINT16;
        else if (dt==4) datatype = V3D_FLOAT32;
        return true;
    }
    else //use Bioformats IO plugin
    {
        QString outfilename;
        if(!call_bioformats_io(imgSrcFile, outfilename))
            return false;

        if (loadImage((char *)qPrintable(outfilename), data1d, sz,  dt))
        {
            if (dt==1) datatype = V3D_UINT8;
            else if (dt==2) datatype = V3D_UINT16;
            else if (dt==4) datatype = V3D_FLOAT32;
            return true;
        }
        else
            return false;
    }

    return false;
}
void DataImporter::cleanData(){
    if (image4d) {delete image4d; image4d = 0;}
}

//double DataImporter::getChannalMinIntensity(V3DLONG channo) //if channo <0 or out of range, then return the in of all channels
//{
//    if (p_vmin && channo>=0 && channo<image4d->getCDim()) return p_vmin[channo];
//    else {V3DLONG tmppos; return minInVector(p_vmin, image4d->getCDim(), tmppos);}
//}

//double DataImporter::getChannalMaxIntensity(V3DLONG channo) //if channo <0 or out of range, then return the max of all channels
//{
//    if (p_vmax && channo>=0 && channo<image4d->getCDim()) return p_vmax[channo];
//    else {V3DLONG tmppos; return maxInVector(p_vmax, image4d->getCDim(), tmppos);}
//}

/***this function should be the one in volimg_proc_declare.h***/
template <class T> bool minMaxInVector(T * p, V3DLONG len, V3DLONG &pos_min, T &minv, V3DLONG &pos_max, T &maxv)
{
    if (!p || len <= 0)
        return false;

    minv = maxv = p[0];
    pos_min = pos_max = 0;
    for (V3DLONG i=1;i<len;i++)
    {
        if (p[i]>maxv)
        {
            maxv = p[i];
            pos_max = i;
        }
        else if (p[i]<minv)
        {
            minv = p[i];
            pos_min = i;
        }
    }

    return true;
}
void DataImporter::updateMinMax(V3DLONG nFrame)
{
    if (image4d)
    {
        //image4d->updateminmaxvalues();

        V3DLONG sx, sy, sz, sc;

        sx = image4d->getXDim();
        sy = image4d->getYDim();
        sz = image4d->getZDim();
        sc = image4d->getCDim();

        if(nFrame<0 || nFrame>=sz) // changed by YuY Feb 8, 2011
            return;

        V3DLONG offsets = nFrame*sx*sy;
        V3DLONG pagesz = sx*sy;

        switch (image4d->getDatatype())
        {
            case V3D_UINT8:
                for(V3DLONG i=0;i<sc;i++)
                {
                    unsigned char minvv,maxvv;
                    V3DLONG tmppos_min, tmppos_max;
                    unsigned char *datahead = (unsigned char *)(image4d->getRawDataAtChannel(i));

                    minMaxInVector(datahead+offsets, pagesz, tmppos_min, minvv, tmppos_max, maxvv);

                    if(p_vmax[i]<maxvv)
                        p_vmax[i] = maxvv;
                    if(p_vmin[i]>minvv)
                        p_vmin[i] = minvv;
                }
                break;

            case V3D_UINT16:
                for(V3DLONG i=0;i<sc;i++)
                {
                    unsigned short int minvv,maxvv;
                    V3DLONG tmppos_min, tmppos_max;
                    unsigned short int *datahead = (unsigned short int *)(image4d->getRawDataAtChannel(i));

                    minMaxInVector(datahead+offsets, pagesz, tmppos_min, minvv, tmppos_max, maxvv);

                    if(p_vmax[i]<maxvv)
                        p_vmax[i] = maxvv;
                    if(p_vmin[i]>minvv)
                        p_vmin[i] = minvv;
                }
                break;

            case V3D_FLOAT32:
                for(V3DLONG i=0;i<sc;i++)
                {
                    float minvv,maxvv;
                    V3DLONG tmppos_min, tmppos_max;
                    float *datahead = (float *)(image4d->getRawDataAtChannel(i));

                    minMaxInVector(datahead+offsets, pagesz, tmppos_min, minvv, tmppos_max, maxvv);

                    if(p_vmax[i]<maxvv)
                        p_vmax[i] = maxvv;
                    if(p_vmin[i]>minvv)
                        p_vmin[i] = minvv;
                }
                break;

            default:
                v3d_msg("Invalid data type found in updateMinMax().");
                return;
        }
    }
}

bool DataImporter::updateminmaxvalues() // update the max min intensity of a single frame (frame 0)
{
    //always delete the two pointers and recreate because if the image is altered in a plugin, the # of color channels may change
    if (p_vmax) {delete []p_vmax; p_vmax=0;}
    if (p_vmin) {delete []p_vmin; p_vmin=0;}

    try
    {
        p_vmax = new double [image4d->getCDim()];
        p_vmin = new double [image4d->getCDim()];
    }
    catch (...)
    {
        v3d_msg("Error happened in allocating memory in updateminmaxvalues().\n");
        if (p_vmax) {delete []p_vmax; p_vmax=0;}
        if (p_vmin) {delete []p_vmin; p_vmin=0;}
        return false;
    }

    V3DLONG i;
    V3DLONG channelPageSize = image4d->getTotalUnitNumberPerChannel();

    switch (image4d->getDatatype())
    {
        case V3D_UINT8:
            for(i=0;i<image4d->getCDim();i++)
            {
                unsigned char minvv,maxvv;
                V3DLONG tmppos_min, tmppos_max;
                unsigned char *datahead = (unsigned char *)image4d->getRawDataAtChannel(i);
                minMaxInVector(datahead, channelPageSize, tmppos_min, minvv, tmppos_max, maxvv);
                p_vmax[i] = maxvv; p_vmin[i] = minvv;
                v3d_msg(QString("channel %1 min=[%2] max=[%3]").arg(i).arg(p_vmin[i]).arg(p_vmax[i]),0);
            }
            break;

        case V3D_UINT16:
            for(i=0;i<image4d->getCDim();i++)
            {
                unsigned short int minvv,maxvv;
                V3DLONG tmppos_min, tmppos_max;
                unsigned short int *datahead = (unsigned short int *)image4d->getRawDataAtChannel(i);
                minMaxInVector(datahead, channelPageSize, tmppos_min, minvv, tmppos_max, maxvv);
                p_vmax[i] = maxvv; p_vmin[i] = minvv;
                v3d_msg(QString("channel %1 min=[%2] max=[%3]").arg(i).arg(p_vmin[i]).arg(p_vmax[i]),0);
            }
            break;

        case V3D_FLOAT32:
            for(i=0;i<image4d->getCDim();i++)
            {
                float minvv,maxvv;

                V3DLONG tmppos_min, tmppos_max;
                float *datahead = (float *)image4d->getRawDataAtChannel(i);

                if (0) //for debugging purpose. 2014-08-22
                {
                    minvv=datahead[0], maxvv=datahead[0];
                    for (V3DLONG myii=1; myii<channelPageSize;myii++)
                    {
                        if (minvv>datahead[myii]) minvv=datahead[myii];
                        else if (maxvv<datahead[myii]) maxvv=datahead[myii];
                    }

                    p_vmax[i] = maxvv; p_vmin[i] = minvv;
                    v3d_msg(QString("channel %1 min=[%2] max=[%3]").arg(i).arg(p_vmin[i]).arg(p_vmax[i]));
                }
                else
                {
                    if (minMaxInVector(datahead, channelPageSize, tmppos_min, minvv, tmppos_max, maxvv))
                    {
                        p_vmax[i] = maxvv; p_vmin[i] = minvv;
                        v3d_msg(QString("channel %1 min=[%2] max=[%3]").arg(i).arg(p_vmin[i]).arg(p_vmax[i]), 0);
                    }
                    else
                    {
                        v3d_msg("fail");
                    }
                }

            }
            break;

        default:
            v3d_msg("Invalid data type found in updateminmaxvalues(). Should never happen, - check with V3D developers.");
            return false;
    }

    return true;
}
