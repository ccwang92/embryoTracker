#include "synquant_simple.h"

//public:
//    Mat zMap, idMap;
//    vector<float> zscore_list;
//    vector<size_t> intensity_levels;
//    float max_exist_zscore;
//    size_t maxZ_intensity_level;
//protected:
//    Mat sortedIdx, parNode;
//    Mat voxSum, voxSumN, areas, areasN, usedN;
//    size_t width, height, zSlice;
//    Mat boxCoor;

using namespace cv;
using namespace std;
//using namespace volproc;
/**
 * @brief synQuantSimple: the major function of synQuant, where the input is the whole FOV data for
 * a thorough cell segmentation
 * @param srcVolume: float intput
 * @param zMap: float output
 * @param idMap: unsigned int output
 * @param p4segVol
 * @param p4odStats
 */
synQuantSimple::synQuantSimple(Mat *_srcVolume, float _src_var, segParameter &p4segVol, odStatsParameter &p4odStats){
    srcVolumeUint8 = _srcVolume;
    imArray = (unsigned char*)srcVolumeUint8->data;
    src_var = _src_var;

    componentTree3d(p4segVol, p4odStats);
    //float * zscore_pointer = & (*zscore_list.begin()); // a pointer to the address of zscore_list.begin()
    zMap = new Mat(srcVolumeUint8->dims, srcVolumeUint8->size, CV_32F, Scalar(0)); //*zscore_list.begin()
    FOREACH_i_ptrMAT(zMap){
        if(outputArray[i] == OBJECT && zscore_list[i] > 1.96){//
            zMap->at<float>(i) = zscore_list[i];
        }
    }
    //ccShowSlice3Dmat(zMap, CV_32F, 3);
    idMap = new Mat(zMap->dims, zMap->size, CV_32S, Scalar(0));
    //Mat b_map = *zMap > 0;
    //ccShowSlice3Dmat(&b_map, CV_8U, 3);
    //cell_num = connectedComponents3d(&b_map, *idMap, 26);
    //ccShowSliceLabelMat(idMap, 3);
    cell_num = floatMap2idMap(zMap, *idMap, 26);
    //ccShowSliceLabelMat(idMap, 0);
//    for(int i = 0; i < idMap->size[2]; i++)
//        ccShowSliceLabelMat(idMap, i);
}
/**
 * @brief synQuantSimple: if the input is a seed for a single cell segmentation
 * @param seed
 */
synQuantSimple::synQuantSimple(singleCellSeed &seed){
    srcVolumeUint8 = &seed.volUint8;
    imArray = (unsigned char*)srcVolumeUint8->data;
    vector<float> var_vals_in_fg = extractValsGivenMask(&seed.varMap, CV_32F, &seed.validSearchAreaMap, 0);
    src_var = (float)vec_mean(var_vals_in_fg);

    idMap = new Mat(srcVolumeUint8->dims, srcVolumeUint8->size, CV_32S, Scalar(0));
    zMap = nullptr; // zmap will not be used if input a seed
    cell_num = 0;
    if(seed.bestFgThreshold > 0){
        maxZ_intensity_level = seed.bestFgThreshold;
    }
    //// the following 5 steps are moved into class cellsegment_main.refineSeed2Region()
//    //// 1. find best threshold to get fgMap and apply some basic morphlogical operation (ordstat4fg.m)
//    // PS: the region from this step will has no overlap with other cell's terriotory
//    // PPS: the cell territory may contain more than one cell, we need verify it using gaps.
//    if(seed.bestFgThreshold < 0){
//        // we have no idea what is the best intensity threshold to get the region
//        // we use this funcion to get this.fgMap
//        cellTerritoryExtractFromSeed(seed, p4odStats);
//    } else{
//        // we have alreay have some idea what is the best intensity level
//        // we use thresholding to get this.fgMap
//        bitwise_and(seed.validSearchAreaMap, seed.volUint8 >= seed.bestFgThreshold, fgMap);
//        // fgMap from simple thresholding may cover >1 cell seeds.
//        removeOtherSeedsInfgMap(seed, p4segVol);
//    }
//    refineCellTerritoryWithSeedRegion(seed, p4segVol);

//    //// 2. update seed's score map based on idMap (fg)
//    normalize(seed.eigMap2d, seed.score2d, 0.001, 1, NORM_MINMAX, CV_32F, fgMap);
//    normalize(seed.eigMap3d, seed.score3d, 0.001, 1, NORM_MINMAX, CV_32F, fgMap);
//    seed.scoreMap = seed.score2d + seed.score3d;
//    //// 3. segment fgMap into idMap based on gaps from principal curvature
//    gapBasedRegionSegment(seed, p4segVol, p4odStats);
//    //// 4. refine the idMap from step 3 based on size and other prior knowledge
//    if(cell_num > 1){
//        bool link_bg2sink = false; // this can force seeds to grow as much as they can
//        Mat input_seeds;
//        idMap->copyTo(input_seeds);
//        //Mat scoreMap = seed.score2d + seed.score3d;
//        regionGrow(&input_seeds, cell_num, *idMap, &seed.scoreMap, &fgMap, //idMap is fully included in fgMap
//                   p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);
//        //TODO: re-assign extra voxels missed in idMap, but contained in this.fgMap
//    }
//    if(cell_num > 0){ // shrink test
//        cellShrinkTest(seed, p4segVol);
//    }
//    removeSmallCC(*idMap, cell_num, p4segVol.min_cell_sz, true);

//    seed.outCell_num = cell_num;
//    idMap->copyTo(seed.outputIdMap); // output
//    seed.bestFgThreshold = (int) maxZ_intensity_level;
//    //// 5. test if the idMap (fg) is too small, if so, enlarge it and re-do the previous steps
//    // PS: this step id done in cellsegment_main.cpp
}


/**
 * @brief cellExtractFromSeed: Currently, we use the exact version of synQuant. But this function should also be able
 * to be implemented using component tree condition on the assumption that the intensity has such trend:
 * inner region > boundary > neighbor. The boundary selection can be fulfilled with the similar way as neighbor selection.
 * @param seed
 * @param p4odStats
 */
void synQuantSimple::cellTerritoryExtractFromSeed(singleCellSeed &seed, odStatsParameter &p4odStats, size_t minSize){
    //Scalar seed_val = mean(seed.volUint8, seed.validSearchAreaMap);// scalar is a vector<double> with length 4 for RGBA data
    //ccShowSlice3Dmat(&seed.volUint8, CV_8U, 3);
    int ub = MAX(30, round(mat_mean(&seed.volUint8, CV_8U, seed.idx_yxz_cropped)));
    int lb = MAX(5, round(mean(seed.volUint8, seed.validSearchAreaMap).val[0]));
    //Mat curr_seed_map = seed.idMap == seed.id; //replaced by seed.seedMap
    if (ub <= lb){
        fgMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
        //qInfo("No valid fg can be found.");
        return;
    }
    vector<float> zscore (ub-lb+1);
    //qInfo("%d thresholds totally: ub:%d, lb:%d.", ub-lb+1, ub, lb);
    Mat otherCellTerritory, valid_cell_territory, cur_reg, cur_valid_nei; // uint8
    int shift_othercell[] = {3,3,1};
    volumeDilate(&seed.otherIdMap, otherCellTerritory, shift_othercell, MORPH_ELLIPSE);
    bitwise_and(otherCellTerritory, seed.seedMap==0, otherCellTerritory); // remove the area covered by seed region
    bitwise_and(otherCellTerritory == 0, seed.validSearchAreaMap, valid_cell_territory); // get the area for current cell and its neighbor
    double seed_sz = (double)seed.idx_yxz.size();
    max_exist_zscore = -1;
    maxZ_intensity_level = -1;
    FOREACH_i(zscore){
        bitwise_and(seed.volUint8 >= (ub-i), valid_cell_territory, cur_reg);
        validSingleRegionExtract(cur_reg, &seed.seedMap, p4odStats.connectInSeedRefine);
        bitwise_and(valid_cell_territory, cur_reg==0, cur_valid_nei);

        size_t reg_sz = overlap_mat_vec(&cur_reg, CV_8U, seed.idx_yxz_cropped, 0);
        if ((reg_sz / seed_sz) < 0.5){
            continue;
        }
        float cur_zscore = debiasedFgBgBandCompare(&cur_reg, &cur_valid_nei, &seed, p4odStats);
        //qInfo("threshold:%d, zscore:%.2f", ub-i, cur_zscore);
        if (max_exist_zscore < cur_zscore){
            maxZ_intensity_level = (ub-i);
            max_exist_zscore = cur_zscore;
        }
    }
    if (max_exist_zscore > 0){
        bitwise_and(seed.volUint8 >= maxZ_intensity_level, valid_cell_territory, fgMap);
        validSingleRegionExtract(fgMap, &seed.seedMap, p4odStats.connectInSeedRefine);
        //qInfo("threshold:%ld, zscore:%.2f", maxZ_intensity_level, max_exist_zscore);
    }else{
        fgMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
        //qInfo("No valid fg can be found.");
        return;
    }
    // test the resultant region
    size_t fg_sz = fgMapSize(&fgMap, CV_8U);
    if (fg_sz < minSize){//isempty_mat_vec(&fgMap, CV_8U, seed.idx_yxz_cropped, 0)){
        fgMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
        //qInfo("No valid fg can be found.");
        return;
    }
}
void synQuantSimple::componentTree3d(segParameter p4segVol, odStatsParameter p4odStats){
    zSlice = srcVolumeUint8->size[2]; // rows x cols x zslices
    width = srcVolumeUint8->size[1];
    height = srcVolumeUint8->size[0];
    size_t nVoxels=width*height*zSlice; // voxels in a 3D image
    size_t nPixels = width*height; //pixels in a single slice
    //imArray = new short[nVoxels];
    //int tmpCnt = 0;

    //    for(int i = 0; i < imArrayIn.length; i++) {
    //        for(int j = 0; j < imArrayIn[0].length; j++)
    //            imArray[tmpCnt++] = imArrayIn[i][j]; // java.util.Arrays
    //    }
    //Create nodes
    long x,y,z, rmder, x0,x2,y0,y2,z0, z2;
    long i,j,k;
    //outputArray = new byte[nVoxels];
    diffN.resize(nVoxels);
    fill(diffN.begin(), diffN.end(), 0);
    voxSumN.resize(nVoxels);
    fill(voxSumN.begin(), voxSumN.end(), 0);
    usedN.resize(nVoxels);
    fill(usedN.begin(), usedN.end(), NOT_USED_AS_N); // or usedNIn[i]
    areas.resize(nVoxels);
    fill(areas.begin(), areas.end(), 1);
    areasN.resize(nVoxels);
    fill(areasN.begin(), areasN.end(), 0);

    voxSum.resize(nVoxels);
    BxCor.resize(nVoxels);//ymin,xmin,zmin,  ymax,xmax,zmax.

    for (i=0; i<nVoxels; i++)
    {
        rmder = i % nPixels;
        z = i / nPixels;
        y=rmder/width;
        x=rmder-y*width;

        voxSum[i] = imArray[i];

        BxCor[i].resize(6);
        BxCor[i][0] = y;
        BxCor[i][1] = x;
        BxCor[i][2] = z;
        BxCor[i][3] = y;
        BxCor[i][4] = x;
        BxCor[i][5] = z;
    }
    parNode.resize(nVoxels);

    //Sort points
    // create a counting array, counts, with a member for
    // each possible discrete value in the input.
    // initialize all counts to 0.
    int maxP = 255; //default for 8-bit image
    int minP = 0;
    size_t nLevels = maxP-minP + 1;
    vector<long> counts(nLevels, 0);
    // for each value in the unsorted array, increment the
    // count in the corresponding element of the count array
    for (i=0; i<nVoxels; i++)
    {
        counts[srcVolumeUint8->at<unsigned char>(i)-minP]++;
    }
    // accumulate the counts - the result is that counts will hold
    // the offset into the sorted array for the value associated with that index
    for (i=1; i<nLevels; i++)
    {
        counts[i] += counts[i-1];
    }
    // store the elements in a new ordered array
    sortedIndex.resize(nVoxels);
    for (i = nVoxels-1; i >= 0; i--)
    {
        // decrementing the counts value ensures duplicate values in A
        // are stored at different indices in sorted.
        sortedIndex[--counts[srcVolumeUint8->at<unsigned char>(i)-minP]] = i;
    }

    //Init nodes
    for (i=0; i<nVoxels; i++)
    {
        parNode[i]=i;
    }
    //Search in decreasing order
    long curNode;
    long adjNode;
    long ii,jj, kk, tmpIdx;
    bool found;
    for (i = nVoxels-1; i >= 0; i--)
    {
        j=sortedIndex[i];
        curNode=j;
        //System.out.println("Image Value"+imArray[j]);
        rmder = j % nPixels;
        z = j / nPixels;
        y=rmder/width;
        x=rmder-y*width;


        found = false;
        /* ROI selection: for the same z-stack, we use 8 neighbors, but for other z-stacks, we only consider two direct neighbors
                * the order of the neighbors are very important
                * we go through the larger neighbors first, then lower ones
                */
        y0=y-1;
        y2=y+1;
        x0=x-1;
        x2=x+1;
        z0=z-1;
        z2=z+1;
        /**for debug*
                if (imArray[j]>0) {
                    System.out.print(j+" " + imArray[j]+"\n");
                }
                if(imArray[j]>=255)//usedN[28+width*28+6*nPixels] != NOT_USED_AS_N)
                    System.out.println(" "+x+" "+y+" "+z);
                */
        //Later neigbours x2,y2
        if(z2<zSlice) {
            k = x+width*y+z2*nPixels;
            if(imArray[k]>=imArray[j])
            {
                adjNode=findNode(k);
                if(curNode!=adjNode)
                {
                    curNode=mergeNodes(adjNode,curNode);
                    found = true;
                }
            }
        }
        if(y2<height)
        {
            k=x+width*y2+z*nPixels;
            if(imArray[k]>=imArray[j])
            {
                adjNode=findNode(k);
                if(curNode!=adjNode)
                {
                    curNode=mergeNodes(adjNode,curNode);
                    found = true;
                }
            }
            if(x2<width)
            {
                k=x2+width*y2+z*nPixels;
                if(imArray[k]>=imArray[j])
                {
                    adjNode=findNode(k);
                    if(curNode!=adjNode)
                    {
                        curNode=mergeNodes(adjNode,curNode);
                        found = true;
                    }
                }
            }
            if(x0>=0)
            {
                k=x0+width*y2+z*nPixels;
                if(imArray[k]>=imArray[j])
                {
                    adjNode=findNode(k);
                    if(curNode!=adjNode)
                    {
                        curNode=mergeNodes(adjNode,curNode);
                        found = true;
                    }

                }
            }
        }
        if(x2<width)
        {
            k=x2+width*y+z*nPixels;
            if(imArray[k]>=imArray[j])
            {
                adjNode=findNode(k);
                if(curNode!=adjNode)
                {
                    curNode=mergeNodes(adjNode,curNode);
                    found = true;
                }
            }
        }
        //Earlier neighbours x0,y0. No need to check =
        if(z0>=0) {
            k = x+width*y+z0*nPixels;
            if(imArray[k]>imArray[j])
            {
                adjNode=findNode(k);
                if(curNode!=adjNode)
                {
                    curNode=mergeNodes(adjNode,curNode);
                    found = true;
                }
            }
        }
        if(x0>=0)
        {
            k=x0+width*y+z*nPixels;
            if(imArray[k]>imArray[j])
            {
                adjNode=findNode(k);
                if(curNode!=adjNode)
                {
                    curNode=mergeNodes(adjNode,curNode);
                    found = true;
                }

            }
        }
        if (y0 >= 0) {
            k = x + width * y0+z*nPixels;
            if (imArray[k] > imArray[j]) {
                adjNode = findNode(k);
                if (curNode != adjNode) {
                    curNode = mergeNodes(adjNode,curNode);
                    found = true;
                }
            }
            if(x2<width)
            {
                k=x2+width*y0+z*nPixels;
                if(imArray[k]>imArray[j])
                {
                    adjNode=findNode(k);
                    if(curNode!=adjNode)
                    {
                        curNode=mergeNodes(adjNode,curNode);
                        found = true;
                    }
                }
            }
            if(x0>=0)
            {
                k=x0+width*y0+z*nPixels;
                if(imArray[k]>imArray[j])
                {
                    adjNode=findNode(k);
                    if(curNode!=adjNode)
                    {
                        curNode=mergeNodes(adjNode,curNode);
                        found = true;
                    }

                }
            }
        }

        if (!found)
        {
            /*****Debug***
                    if(j==13979) {
                        System.out.print("neighbor: "+voxSumN[j]+" "+areasN[j]+" self: "+voxSum[j]+" "+areas[j]+" "+"\n");
                    }***/
            y0= MAX(y-1,0);
            y2= MIN(y+1, height-1);
            x0= MAX(x-1,0);
            x2= MIN(x+1,width-1);
            z0= MAX(z-1,0);
            z2= MIN(z+1,zSlice-1);
            // for neighboring pixels' value we consider 26 neighbors
            for (ii=z2;ii>=z0;ii--) {
                for (jj=y2;jj>=y0;jj--) {
                    for(kk=x2;kk>=x0;kk--) {
                        if( ii==z & jj==y & kk==x)
                            continue;
                        tmpIdx = kk+width*jj+ii*nPixels;
                        if (usedN[tmpIdx] == NOT_USED_AS_N)
                        {
                            voxSumN[j] += imArray[tmpIdx];
                            areasN[j]++;
                            usedN[tmpIdx] = USED_AS_N_ONCE;
                            /*if(j==13979) {
                                            System.out.print("neighbor val: "+imArray[tmpIdx]+" "+ii +" "+jj+" "+kk+"\n");
                                        }*/
                        }
                    }
                }
            }
            usedN[j] = USED_AS_N_MORE;
            diffN[j] = voxSum[j]/(double)areas[j];
            if (areasN[j] > 0)
                diffN[j] -= voxSumN[j]/(double)areasN[j];
        }

    }
    outputArray.resize(nVoxels); // label the output
    fill(outputArray.begin(), outputArray.end(), UNDEFINED);
    zscore_list.resize(nVoxels);
    for (i=0; i<nVoxels; i++)
    {
        rmder = i % nPixels;
        z = i / nPixels;
        y=rmder/width;
        x=rmder-y*width;
        double LH = (double)BxCor[i][4]-BxCor[i][1]+1;
        double LW = (double)BxCor[i][3]-BxCor[i][0]+1;
        double LZ = (double)BxCor[i][5]-BxCor[i][2]+1;
        double ratio = LH>LW? LH/LW: LW/LH;

        if(areas[i]>=p4segVol.max_cell_sz){
            zscore_list[i] = -1;
            outputArray[i] = NOT_OBJECT; // no need to update the object label
        }
        if(areas[i]<p4segVol.min_cell_sz || ratio>p4segVol.max_WHRatio || (areas[i]/(double)(LH*LW*LZ))<p4segVol.min_fill){
            zscore_list[i] = -1;
        }else{
            zscore_list[i] = zscoreCal(diffN[i], areas[i], areasN[i]);
            valid_zscore.push_back(zscore_list[i]);
            valid_zscore_idx.push_back(i);
        }
    }
//    Mat tmp = Mat(3, srcVolumeUint8->size, CV_32F, zscore_list.data());
//    ccShowSlice3Dmat(tmp, CV_32F);
    // label the objects
    objLabel_descending();
}

/*Label object or not: here we simply  test the size  to decide object or not*/
void synQuantSimple::objLabel(size_t minSize ,size_t maxSize)
{
    long i,j;
    size_t nVoxels = sortedIndex.size();
    for (i = nVoxels-1; i >= 0; i--)
    {
        j=sortedIndex[i];
        if (areas[j]<=maxSize && areas[j]>=minSize)
            outputArray[j] = OBJECT;
        else
            outputArray[j] = NOT_OBJECT;
    }
}
size_t synQuantSimple::findNode(size_t e)
{
    if(parNode[e]!=e)
    {
        size_t root = findNode(parNode[e]);
        //parNode[e] = root; //This cannot be used here
        return root;
    }
    else
    {
        return e;
    }
}
size_t synQuantSimple::mergeNodes(size_t e1,size_t e2)/*e1 adjacent node, e2 current node*/
{
    //e1-adjacent; e2-current
    size_t res;
    size_t m;

    if(imArray[e1]==imArray[e2])
    {
        res=max(e1,e2);
        m=min(e1,e2);
    }
    else
    {
        res=e2;
        m=e1;
    }
    size_t curNeiCnt = areasN[res];
    if (curNeiCnt==53)
        curNeiCnt = 53;
    /*****Debug****
        if(res==13979) {
            System.out.print("neighbor: "+voxSumN[res]+" "+areasN[res]+" self: "+voxSum[res]+" "+areas[res]+" "+"\n");
        }*/
    //Compute new neighbours
    long z = e2 / (width*height);
    long rmder = e2 % (width*height);
    long y=rmder/width;
    long x=rmder-y*width;

    long y0=MAX(y-1,0);
    long y2=MIN(y+1, height-1);
    long x0=MAX(x-1,0);
    long x2=MIN(x+1,width-1);
    long z0=MAX(z-1,0);
    long z2=MIN(z+1,zSlice-1);
    // for neighboring pixels' value we consider 26 neighbors
    //System.out.print("Before Merging"+voxSumN[res]+" "+areasN[res]+" "+"\n");
    long ii, jj, kk, tmpIdx;
    for (ii=z2;ii>=z0;ii--) {
        for (jj=y2;jj>=y0;jj--) {
            for(kk=x2;kk>=x0;kk--) {
                if( ii==z && jj==y && kk==x){
                    continue;
                }
                tmpIdx = kk+width*jj+ii*(width*height);
                if (usedN[tmpIdx] == NOT_USED_AS_N)
                {
                    voxSumN[res] += imArray[tmpIdx];
                    areasN[res]++;
                    usedN[tmpIdx] = USED_AS_N_ONCE;
                }
            }
        }
    }
    /*****Debug***
        if(res==136)
            System.out.print("e2: "+areasN[e2]+ "e1: "+areasN[e1]+"\n");
        **/
    areasN[res] += areasN[m];
    voxSumN[res] += voxSumN[m];
    if (usedN[e2] == USED_AS_N_ONCE) // e2 ever be used as neighbors, now we need to remove them
    {
        areasN[res] -= areas[e2];
        voxSumN[res] -= voxSum[e2];
    }
    /*****Debug***
        if(res==136)
            System.out.print("Before Merging res: "+curNeiCnt+", but after Merging res: "+areasN[res]+"\n");
        */
    areas[res] += areas[m];
    voxSum[res] += voxSum[m];
    parNode[m]=res;

    usedN[e2] = USED_AS_N_MORE;

    diffN[res] = voxSum[res]/(double)areas[res];

    if (areasN[res] > 0)
        diffN[res] -= voxSumN[res]/(double)areasN[res];

    //System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
    //System.out.println("BC:" +BxCor[m][3]+" "+BxCor[m][2]+" "+BxCor[m][1]+" "+BxCor[m][0]);
    if (BxCor[res][0] > BxCor[m][0])
        BxCor[res][0] = BxCor[m][0];
    if (BxCor[res][1] > BxCor[m][1])
        BxCor[res][1] = BxCor[m][1];
    if (BxCor[res][2] > BxCor[m][2])
        BxCor[res][2] = BxCor[m][2];
    if (BxCor[res][3] < BxCor[m][3])
        BxCor[res][3] = BxCor[m][3];
    if (BxCor[res][4] < BxCor[m][4])
        BxCor[res][4] = BxCor[m][4];
    if (BxCor[res][5] < BxCor[m][5])
        BxCor[res][5] = BxCor[m][5];
    //System.out.println("BC:" +BxCor[res][3]+" "+BxCor[res][2]+" "+BxCor[res][1]+" "+BxCor[res][0]);
    return res;
}
/*Calculate the z-score of each pixel*/
float synQuantSimple::zscoreCal(float t0, size_t M/*in*/, size_t N/*nei*/){
    float mu, sigma;
    nonOV_truncatedGauss(M, N, mu, sigma);
    mu = mu*sqrt(src_var);
    sigma = sigma*sqrt(src_var);
    t0 = t0/sqrt(src_var);
    return (t0-mu)/sigma;
}

/**
 * @brief objLabel_zscore: label the objects by z-score threshold set by users
 * @param zscore_thres
 */
void synQuantSimple::objLabel_zscore(float zscore_thres) {
//    int i,j;
//    for (i=nVoxels-1; i>=0; i--)
//    {
//        j=sortedIndex[i];
//        if (outputArray[j] == UNDEFINED)
//        {
//            int e = j;
//            while(zscore[e] < zscore_thres && outputArray[e] == UNDEFINED) {
//                e = parNode[e];
//            }
//            if (zscore[e] >= zscore_thres) { // this region is good
//                double cur_zscore = zscore[e];
//                double cur_label = outputArray[e];
//                while(imArray[e] == imArray[parNode[e]] & outputArray[e] == UNDEFINED)
//                    e = parNode[e];
//                if (cur_label == UNDEFINED) { // this region haven't been labeled
//                    int e1 = j;
//                    while(e1 != e)
//                    {
//                        outputArray[e1] = OBJECT;
//                        zscore[e1] = cur_zscore;
//                        e1 = parNode[e1];
//                    }
//                    if (outputArray[e1] == UNDEFINED) {
//                        outputArray[e1] = OBJECT;
//                        zscore[e1] = cur_zscore;
//                    }
//                    e1 = parNode[e];
//                    while(e1 != parNode[e1] & outputArray[e1] == UNDEFINED)
//                    {
//                        outputArray[e1] = NOT_OBJECT;
//                        zscore[e1] = -1;
//                        e1 = parNode[e1];
//                    }
//                    if (outputArray[e1] == UNDEFINED) {
//                        outputArray[e1] = NOT_OBJECT;
//                        zscore[e1] = -1;
//                    }
//                    if (outputArray[e1] == OBJECT) { // we did wrong thing if the root is OBJECT
//                        // under such condition, we need to correct previous from j to e1
//                        int e2 = j;
//                        while(e2 != e1)
//                        {
//                            outputArray[e2] = OBJECT;
//                            zscore[e2] = zscore[e1];
//                            e2 = parNode[e2];
//                        }
//                    }
//                }
//                else{ // this region has been labeled (should only have NOT_OBJECT label)
//                    int e1 = j;
//                    while(e1 != e) // label j to e as the same
//                    {
//                        outputArray[e1] = cur_label;
//                        zscore[e1] = cur_zscore;
//                        e1 = parNode[e1];
//                    }
//                }
//            }else { // this region is bad
//                int e1 = j;
//                while(e1 != e) // label j to e as the same
//                {
//                    outputArray[e1] = outputArray[e];
//                    zscore[e1] = zscore[e];
//                    e1 = parNode[e1];
//                }
//            }
//        }
//    }
}

/**
 * @brief synQuantSimple::objLabel_fdr: label the objects by finding the
 * best level based on z-score
 */
void synQuantSimple::processVoxLabel(size_t j){
    if (outputArray[j] == UNDEFINED)
    {
        size_t e = j;
        //System.out.println("Idx "+j+" Image Value"+imArray[j]+" "+zscore[j]);
        while(imArray[e] == imArray[parNode[e]] && outputArray[e] == UNDEFINED)
            e = parNode[e];
        if (outputArray[e] == UNDEFINED)
        {
            size_t e1 = parNode[e];
            //					if (e1!=e && zscore[e]>0) {
            //						outputArray[e] = OBJECT;
            //						continue;
            //					}
            // find if parents contains a valid OBEJCT or NOT_OBJECT label
            byte ObjOrNot = NOT_OBJECT;
            double tmpZscore = -1;
            while(e1 != parNode[e1])
            {
                if (outputArray[e1]!=UNDEFINED) {
                    ObjOrNot = outputArray[e1];
                    tmpZscore = zscore_list[e1];
                    break;
                }
                e1 = parNode[e1];
            }


            if (ObjOrNot==OBJECT) { //if there is a OBEJCT label, assign j to e1 to object
                e1 = j;
                while(e1 != parNode[e1] && outputArray[e1] == UNDEFINED){
                    outputArray[e1] = OBJECT;
                    //diffN[e1] = diffN[e];
                    zscore_list[e1] = tmpZscore;
                    e1 = parNode[e1];
                }
                if (outputArray[e1] == UNDEFINED) { // should not use; in case loop till the last pixel
                    outputArray[e1] = OBJECT;
                    zscore_list[e1] = tmpZscore;
                }
                //						if (outputArray[e] == UNDEFINED) {
                //							outputArray[e] = OBJECT;
                //							zscore[e] = tmpZscore;
                //						}
            }
            else { // label parNode[e] to e1 NOT_OBJECT, label j to e object
                e1 = parNode[e];
                while(e1 != parNode[e1] && outputArray[e1] == UNDEFINED){//may lost one pixel; do-while may better
                    outputArray[e1] = NOT_OBJECT;
                    //diffN[e1] = diffN[e];
                    zscore_list[e1] = -1;
                    e1 = parNode[e1];
                }
                if (outputArray[e1] == UNDEFINED) {
                    outputArray[e1] = NOT_OBJECT;
                    zscore_list[e1] = -1;
                }
                // only e and j is enough
                outputArray[e] = OBJECT;
                zscore_list[e] = zscore_list[j];
                outputArray[j] = OBJECT;
                zscore_list[j] = zscore_list[j];
            }

            //findBestLevel(j, -1);
        }
        else
        {
            size_t e1 = j;
            while(e1 != e)
            {
                outputArray[e1] = outputArray[e];
                //diffN[e1] = diffN[e];
                zscore_list[e1] = zscore_list[e];
                e1 = parNode[e1];
            }
        }
    }
}
void synQuantSimple::objLabel_descending() {
    long j;
    // sort all pixels with their zscore values with descending order
    vector<size_t> zScoreSortedIdx = sort_indexes(valid_zscore, false, 0); //false:descend, 0:start from 0
    // process the voxel with valid zscore first
    FOREACH_i(zScoreSortedIdx)
    {
        j=valid_zscore_idx[zScoreSortedIdx[i]];
        processVoxLabel(j);
    }
    // process other voxels, no need to avoid voxels processed in valid_zscore_idx.
    FOREACH_i(zscore_list)
    {
        processVoxLabel(i);
    }
}
/**
 * @brief synQuantSimple::debiasedFgBgBandCompare
 * @param cur_reg
 * @param seed
 * @param p4odStats
 * @param debiasMethod
 * @return
 */
float synQuantSimple::debiasedFgBgBandCompare(Mat *cur_reg, Mat *validNei, singleCellSeed *seed,
                                                    odStatsParameter p4odStats){
    Mat cur_reg_dilate, cur_reg_erode;
    int dilate_size[] = {p4odStats.roundNum4fgbgCompare, p4odStats.roundNum4fgbgCompare, 0};
    volumeDilate(cur_reg, cur_reg_dilate, dilate_size, MORPH_ELLIPSE);
    //volumeErode(cur_reg, cur_reg_erode, dilate_size, MORPH_ELLIPSE);
    //ccShowSlice3Dmat(&cur_reg_dilate, CV_8U);
    Mat fg_center, fg_band, bg_band;
    //bitwise_and(*cur_reg, cur_reg_erode == 0, fg_band);
    bitwise_and(*validNei, cur_reg_dilate, bg_band);
    //ccShowSlice3Dmat(&bg_band, CV_8U);
    volumeDilate(&bg_band, fg_band, dilate_size, MORPH_ELLIPSE);
    //ccShowSlice3Dmat(&fg_band, CV_8U);
    bitwise_and(fg_band, *cur_reg, fg_band);
    //ccShowSlice3Dmat(&fg_band, CV_8U);
    fg_center = seed->validSearchAreaMap - fg_band - bg_band; // key here
    vector<float> fg_band_vals, bg_band_vals, fg_center_vals;
    fg_band_vals = extractValsGivenMask(&seed->volUint8, CV_8U, &fg_band, 0);
    bg_band_vals = extractValsGivenMask(&seed->volUint8, CV_8U, &bg_band, 0);
    fg_center_vals = extractValsGivenMask(&seed->volUint8, CV_8U, &fg_center, 0);

    float mu, sigma, zscore = 0;
    if (p4odStats.fgSignificanceTestWay == KSEC){
        //fg_center_vals.resize(3000);
        orderStatsKSection(fg_band_vals, bg_band_vals, fg_center_vals, mu, sigma);
        //qInfo("mu:%.4f, sigma:%.4f", mu, sigma);
        float sum_stats = vec_mean(fg_band_vals) - vec_mean(bg_band_vals);
        zscore = (sum_stats - mu*sqrt(src_var))/(sigma*sqrt(src_var));
    }else{
        //TODO: other ways haven't been implemented here
    }
    return zscore;
}
