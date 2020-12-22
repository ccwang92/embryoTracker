#ifndef CELLSEGMENT_MAIN_H
#define CELLSEGMENT_MAIN_H
#include <opencv2/core.hpp> //basic building blocks of opencv
#include <opencv2/imgcodecs.hpp> // image io
#include <opencv2/highgui.hpp> //image display
//#define VOLUME_WH_MAX 100000
using namespace cv;
using namespace std;
struct segParameter {
    union
    {
        float min_intensity;
        long min_cell_sz;
        long max_cell_sz;
        float min_fill;
        long max_WHRatio;
        float noise_estimate_ratio;
        float fdr;
        float min_zscore;
    };
};

struct odStatsParameter {
    union
    {
        int gap4fgbgCompare; // default is 0
    };
};

class cellSegmentMain
{
public:
    cellSegmentMain(unsigned char *data_grayim4d, int _data_type, long buffSize[5]/*(x,y,z,c,t)*/);
    void cellSegmentSingleFrame(unsigned char *data_grayim3d, size_t curr_frame);

protected:
    string debug_folder;
    string default_name;
    int data_type;
    int * data_rows_cols_slices;
    long time_points = 0;
    //float min_intensity = 0.0;
    segParameter p4segVol;
    odStatsParameter p4odStats;
    vector<Mat> cell_label_maps;
    vector<Mat> principalCurv2d;
    vector<Mat> principalCurv3d;
    vector<Mat> varMaps;
    vector<size_t> number_cells;
    vector<vector<size_t>> voxIdxList; // cell voxIdx list

    void init_parameter(){
        string debug_folder = "/home/ccw/Desktop/embryo_res_folder/";
        string default_name = debug_folder + "test.tiff";
        p4segVol.min_intensity = 0.0;
        p4segVol.fdr = .05;
        p4segVol.min_cell_sz =100;
        p4segVol.max_cell_sz = 3000;
        p4segVol.min_fill = 0.0001;
        p4segVol.max_WHRatio = 100;

        p4odStats.gap4fgbgCompare = 0;
    }
};

#endif // CELLSEGMENT_MAIN_H
