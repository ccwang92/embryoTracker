#include "data_importer.h"
#include <stdio.h>
#include "dialog_import_im_sequence.h"

#include "src_3rd/basic_c_fun/stackutil.h"
//#include "src_3rd/basic_c_fun/volimg_proc_declare.h"
#include "src_3rd/io/io_bioformats.h"
#include "src_3rd/basic_c_fun/v3d_message.h"
//#include "v3d_message.h"

/**
** brief: main function of loading tiff image sequence
*/
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
            //v3d_msg("File loaded!");
            return true;
        }
        else
            v3d_msg("Fail loading!");
    }
    else
        v3d_msg("Error file name!");
    return false;
}
/**
** brief: load tiff images in the same folder
*/
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

/**
** brief: tranverse images in current folder
*/
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

    if (debugMode){
        qDebug("Currently we are debuging the segment/track algorithm, no pop-up dialog\n");
        packtype = TIME_PACK_C;
        return myList;
    }
    //Get the tiff image sequence by usr interaction

    ImportImgPara p;
    p.countImg = myList.size();

    //import_images_tool_Dialog d(curFilePath);
    dialogImportImSequence d;
    //need to update the information based on the current myList info. 131017
    d.num_im_SpBox->setMaximum(p.countImg);
    d.num_im_SpBox->setValue(p.countImg);
    d.num_im_SpBox->setMinimum(p.countImg);

    d.start_id_SpBox->setMaximum(p.countImg);
    d.start_id_SpBox->setValue(1);

    d.increment_SpBox->setMaximum(p.countImg);

    d.end_id_SpBox->setMaximum(p.countImg);
    d.end_id_SpBox->setValue(p.countImg);

    int res = d.exec();
    if (res==QDialog::Accepted)
    {
        d.fetchData(&p);// get the images for loading

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
/**
** brief: load single tiff image
*/
bool DataImporter::readSingleImageFile(char *imgSrcFile, unsigned char * & data1d, V3DLONG * & sz, ImagePixelType & datatype)
{
    datatype = V3D_UNKNOWN;
    int dt = 0;
    if (loadImage(imgSrcFile, data1d, sz,  dt))
    {
        if (dt==1) datatype = V3D_UINT8;
        else if (dt==2) datatype = V3D_UINT16;
        else if (dt==4) datatype = V3D_FLOAT32;


        if(!image4d->getDatatype())
        {
            image4d->setDatatype(datatype);
            qDebug("Data type is %d\n", image4d->getDatatype());
        }

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

            if(!image4d->getDatatype()) {
                image4d->setDatatype(datatype);
                qDebug("Data type is %d\n", image4d->getDatatype());
            }

            return true;
        }
        else
            return false;
    }

    return false;
}
/**
** brief: destructor
*/
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
        p_vmax = new double [image4d->getCDim()]; // at first the data was ordered by channel
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

/**
 * /brief, update the input data to color image; or the counter-transfer
 **/
void DataImporter::rgba3d_r2gray(RGBA8* rgbaBuf, V3DLONG bufSize[5])
{
    if (rgbaBuf==0 || bufSize==0)
        return;

    // copy R to G,B
    V3DLONG imageX, imageY, imageZ, imageC, imageT;
    {
        imageX = bufSize[0];
        imageY = bufSize[1];
        imageZ = bufSize[2];
        imageC = 4;
        imageT = bufSize[4];
    }
    if (imageX*imageY*imageZ==0)
        return;

    V3DLONG ot;
    V3DLONG ox, oy, oz;
    for (ot=0; ot<imageT; ot++)
    for (oz = 0; oz < imageZ; oz++)
    for (oy = 0; oy < imageY; oy++)
    for (ox = 0; ox < imageX; ox++)
        {
            RGBA8 rgba;
            V3DLONG p = ot*(imageZ*imageY*imageX) + oz*(imageY*imageX) + oy*(imageX) + ox;

            rgba = rgbaBuf[p];
            rgbaBuf[p].g = rgba.r;
            rgbaBuf[p].b = rgba.r;
        }
}

// 090602 RZC: inline cannot be used in *.cpp but in *.h
// Please see the below comment for 'sampling3dUNIT8_2' in relation to this function.
float DataImporter::sampling3dUINT8(Image4DProxy<Image4DSimple>& img4dp,
        V3DLONG c,
        V3DLONG x, V3DLONG y, V3DLONG z, V3DLONG dx, V3DLONG dy, V3DLONG dz)
{
    V3DLONG dim1=img4dp.sx; V3DLONG dim2=img4dp.sy; V3DLONG dim3=img4dp.sz;
    float avg = 0;
    float d = (dx*dy*dz);
    if (d>0 && x>=0 && y>=0 && z>=0 && x+dx<=dim1 && y+dy<=dim2 && z+dz<=dim3)
    {
        //float s = 0;
        V3DLONG xi,yi,zi;
        for (zi=0; zi<dz; zi++)
        for (yi=0; yi<dy; yi++)
        for (xi=0; xi<dx; xi++)
        {
            //float w = MAX( (2-ABS(xi-0.5*dx)*2.0/dx), MAX( (2-ABS(yi-0.5*dy)*2.0/dy), (2-ABS(zi-0.5*dz)*2.0/dz) )); //090712
            //s += w;
            avg += img4dp.value8bit_at(x,y,z,c);// *w;
        }
        //d = s;
        avg /= d;
    }
    return avg;
}

// There is a mystery as to why the above function, sampling3dUNIT8, uses a hard-coded 'img4dp.value8bit_at(x,y,z,c)' rather than
// sampling using the zi,yi,xi indices. Since it does not appear to use the indices, it is equivalent to the below implementation,
// which speeds-up the parent function quite a bit. Until we can resolve what is going on, we are temporarily using this
// fast version.
float DataImporter::sampling3dUINT8_2(Image4DProxy<Image4DSimple>& img4dp,
                V3DLONG c,
                V3DLONG x, V3DLONG y, V3DLONG z, V3DLONG dx, V3DLONG dy, V3DLONG dz, V3DLONG dxyz)
{
        V3DLONG dim1=img4dp.sx; V3DLONG dim2=img4dp.sy; V3DLONG dim3=img4dp.sz;
        float avg = 0;
        if (dxyz>0 && x>=0 && y>=0 && z>=0 && x+dx<=dim1 && y+dy<=dim2 && z+dz<=dim3)
        {
                avg = img4dp.value8bit_at(x,y,z,c);
        }
        return avg;
}

void DataImporter::data4dp_to_rgba3d(Image4DProxy<Image4DSimple>& img4dp, V3DLONG dim5,
        V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
        V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
        RGBA8* rgbaBuf, V3DLONG bufSize[5])
{
    if (rgbaBuf==0 || bufSize==0)
        return;

    //there may be a memory issue? by PHC 20110122
    //	if (img4dp.su!=1)
//	{
//		v3d_msg("Your data is not 8bit. Now this data4dp_to_rgba3d(0 function supports only 8bit data.");
//		return;
//	}

    V3DLONG dim1=img4dp.sx; V3DLONG dim2=img4dp.sy; V3DLONG dim3=img4dp.sz;
    V3DLONG dim4=img4dp.sc;
    #define SAMPLE(it, ic, ix,iy,iz, dx,dy,dz) \
                (unsigned char)sampling3dUINT8( img4dp, (it*dim4/imageT + ic), \
                                                ix, iy, iz, dx, dy, dz )

    #define SAMPLE2(si, ix,iy,iz, dx,dy,dz, dxyz) \
        (unsigned char)sampling3dUINT8_2( img4dp, si, ix, iy, iz, dx, dy, dz, dxyz)

    // only convert 1<=dim4<=4 ==> RGBA
    V3DLONG imageX, imageY, imageZ, imageC, imageT;
    {
        imageX = bufSize[0];
        imageY = bufSize[1];
        imageZ = bufSize[2];
                imageC = MIN(4, size4); // <=4
        imageT = bufSize[4];
    }
    if (imageX*imageY*imageZ*imageC*imageT==0)
        return;

    float sx, sy, sz;
    V3DLONG dx, dy, dz;
    sx = float(size1) / imageX;
    sy = float(size2) / imageY;
    sz = float(size3) / imageZ;
    dx = V3DLONG(sx);
    dy = V3DLONG(sy);
    dz = V3DLONG(sz);
        V3DLONG dxyz = dx*dy*dz;
    MESSAGE_ASSERT(dx*dy*dz >=1); //down sampling

    V3DLONG ot;
    V3DLONG ox, oy, oz;
    V3DLONG ix, iy, iz;

        V3DLONG otOffset, ozOffset, oyOffset, oxOffset;

        V3DLONG SAM0, SAM1, SAM2, SAM3;


        for (ot=0; ot<imageT; ot++) {
            SAM0 = ot*dim4/imageT + 0;
            SAM1 = SAM0+1;
            SAM2 = SAM0+2;
            SAM3 = SAM0+3;
            otOffset=ot*(imageZ*imageY*imageX);
            for (oz = 0; oz < imageZ; oz++) {
                ozOffset=oz*(imageY*imageX);
                iz = start3+ CLAMP(0,dim3-1, IROUND(oz*sz));
                for (oy = 0; oy < imageY; oy++) {
                    oyOffset=oy*imageX;
                    oxOffset=otOffset+ozOffset+oyOffset;
                    iy = start2+ CLAMP(0,dim2-1, IROUND(oy*sy));

                    if (imageC==1) {
                        for (ox=0;ox<imageX;ox++) {
                            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
                            RGBA8 rgba;
                            rgba.r = SAMPLE2(SAM0, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.g = 0;
                            rgba.b = 0;
                            float t = (0.f + rgba.r + rgba.g + rgba.b);
                            rgba.a = (unsigned char)t;
                            rgbaBuf[oxOffset++] = rgba;
                        }
                    }

                    if (imageC==2) {
                        for (ox=0;ox<imageX;ox++) {
                            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
                            RGBA8 rgba;
                            rgba.r = SAMPLE2(SAM0, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.g = SAMPLE2(SAM1, ix,iy,iz, dx,dy,dz, dxyz);;
                            rgba.b = 0;
                            float t = (0.f + rgba.r + rgba.g + rgba.b)/2.0;
                            rgba.a = (unsigned char)t;
                            rgbaBuf[oxOffset++] = rgba;
                        }
                    }

                    if (imageC==3) {
                        for (ox=0;ox<imageX;ox++) {
                            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
                            RGBA8 rgba;
                            rgba.r = SAMPLE2(SAM0, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.g = SAMPLE2(SAM1, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.b = SAMPLE2(SAM2, ix,iy,iz, dx,dy,dz, dxyz);
                            float t = (0.f + rgba.r + rgba.g + rgba.b)/3.0;
                            rgba.a = (unsigned char)t;
                            rgbaBuf[oxOffset++] = rgba;
                        }
                    }

                    if (imageC>=4) {
                        for (ox=0;ox<imageX;ox++) {
                            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
                            RGBA8 rgba;
                            rgba.r = SAMPLE2(SAM0, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.g = SAMPLE2(SAM1, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.b = SAMPLE2(SAM2, ix,iy,iz, dx,dy,dz, dxyz);
                            rgba.a = SAMPLE2(SAM3, ix,iy,iz, dx,dy,dz, dxyz);
                            rgbaBuf[oxOffset++] = rgba;
                        }
                    }
                }
            }
        }
}

//void DataImporter::data4dp_to_rgba3d(unsigned char* data4dp, V3DLONG dim1, V3DLONG dim2, V3DLONG dim3, V3DLONG dim4, V3DLONG dim5,
//        V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
//        V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
//        RGBA8* rgbaBuf, V3DLONG bufSize[5])
//{
//    if (data4dp==0 || rgbaBuf==0 || bufSize==0)
//        return;

//    //V3DLONG block_size = dim3*dim2*dim1;
//    #define SAMPLE(it, ic, ix,iy,iz, dx,dy,dz) \
//                (unsigned char)sampling3dAllTypes( data4dp+ (it*dim4 + ic)*(dim3*dim2*dim1), \
//                                                dim1, dim2, dim3, ix, iy, iz, dx, dy, dz )

//    // only convert 1<=dim4<=4 ==> RGBA
//    V3DLONG imageX, imageY, imageZ, imageC, imageT;
//    {
//        imageX = bufSize[0];
//        imageY = bufSize[1];
//        imageZ = bufSize[2];
//        imageC = MIN(4, size4); // <=4
//        imageT = bufSize[4];
//    }
//    if (imageX*imageY*imageZ*imageC*imageT==0)
//        return;

//    float sx, sy, sz;
//    V3DLONG dx, dy, dz;
//    sx = float(size1) / imageX;
//    sy = float(size2) / imageY;
//    sz = float(size3) / imageZ;
//    dx = V3DLONG(sx);
//    dy = V3DLONG(sy);
//    dz = V3DLONG(sz);
//    MESSAGE_ASSERT(dx*dy*dz >=1); //down sampling

//    V3DLONG ot;
//    V3DLONG ox, oy, oz;
//    V3DLONG ix, iy, iz;
//    for (ot=0; ot<imageT; ot++)
//    for (oz = 0; oz < imageZ; oz++)
//    for (oy = 0; oy < imageY; oy++)
//    for (ox = 0; ox < imageX; ox++)
//        {
//            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
//            iy = start2+ CLAMP(0,dim2-1, IROUND(oy*sy));
//            iz = start3+ CLAMP(0,dim3-1, IROUND(oz*sz));

//            RGBA8 rgba;

//            if (imageC >= 1) {
//                rgba.r = SAMPLE(ot, 0, ix,iy,iz, dx,dy,dz);
//            } else {
//                rgba.r = 0;
//            }

//            if (imageC >= 2) {
//                rgba.g = SAMPLE(ot, 1, ix,iy,iz, dx,dy,dz);
//            } else {
//                rgba.g = 0;
//            }

//            if (imageC >= 3) {
//                rgba.b = SAMPLE(ot, 2, ix,iy,iz, dx,dy,dz);
//            } else {
//                rgba.b = 0;
//            }

//            if (imageC >= 4) {
//                rgba.a = SAMPLE(ot, 3, ix,iy,iz, dx,dy,dz);
//            } else {
//                float t = //MAX(rgba.r, MAX(rgba.g, rgba.b));
//                            ((0.f + rgba.r + rgba.g + rgba.b) / imageC);
//                rgba.a = (unsigned char)t;
//                            //(unsigned char)(t*t/255);
//                            //(unsigned char)(sqrt(t/255)*255);
//            }

//            rgbaBuf[ot*(imageZ*imageY*imageX) + oz*(imageY*imageX) + oy*(imageX) + ox] = rgba;
//        }
//}
