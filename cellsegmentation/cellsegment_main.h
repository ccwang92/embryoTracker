#ifndef CELLSEGMENT_MAIN_H
#define CELLSEGMENT_MAIN_H

//#define VOLUME_WH_MAX 100000

//enum data_type {UINT8 = 1, UINT16 = 2, UNKONWN = 4};
class cellSegmentMain
{
public:
    cellSegmentMain(unsigned char *data_grayim4d, int _data_type, long buffSize[5]/*(x,y,z,c,t)*/);
    void cellSegmentSingleFrame(unsigned char *data_grayim3d);

protected:
    int data_type;
    int * data_rows_cols_slices;
    long time_points = 0;
};

#endif // CELLSEGMENT_MAIN_H
