#include "synquant_simple.h"
#include "img_basic_proc.h"
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

/**
 * @brief synQuantSimple: the major function of synQuant
 * @param srcVolume: float intput
 * @param zMap: float output
 * @param idMap: unsigned int output
 * @param p4segVol
 * @param p4odStats
 */
synQuantSimple::synQuantSimple(Mat *_srcVolume, float _src_var, segParameter p4segVol, odStatsParameter p4odStats){
    srcVolumeUint8 = _srcVolume;
    imArray = (unsigned char*)srcVolumeUint8->data;
    src_var = _src_var;

    componentTree3d(p4segVol, p4odStats);
    float * zscore_pointer = & (*zscore_list.begin()); // a pointer to the address of zscore_list.begin()
    zMap = new Mat(srcVolumeUint8->dims, srcVolumeUint8->size, CV_32F, zscore_pointer); //*zscore_list.begin()
    floatMap2idMap(zMap, *idMap, 26);
}

void synQuantSimple::componentTree3d(segParameter p4segVol, odStatsParameter p4odStats){
    zSlice = srcVolumeUint8->size[2]; // rows x cols x zslices
    width = srcVolumeUint8->size[2];
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
    size_t i,j,k;
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
    size_t counts[nLevels];
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
    size_t curNode;
    size_t adjNode;
    size_t ii,jj, kk, tmpIdx;
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
    // label the object
    objLabel_descending();
}

/*Label object or not: here we simply  test the size  to decide object or not*/
void synQuantSimple::objLabel(size_t minSize ,size_t maxSize)
{
    size_t i,j;
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
    size_t z = e2 / (width*height);
    size_t rmder = e2 % (width*height);
    size_t y=rmder/width;
    size_t x=rmder-y*width;

    size_t y0=max(y-1,(size_t)0);
    size_t y2=min(y+1, height-1);
    size_t x0=max(x-1,(size_t)0);
    size_t x2=min(x+1,width-1);
    size_t z0=max(z-1,(size_t)0);
    size_t z2=min(z+1,zSlice-1);
    // for neighboring pixels' value we consider 26 neighbors
    //System.out.print("Before Merging"+voxSumN[res]+" "+areasN[res]+" "+"\n");
    size_t ii, jj, kk, tmpIdx;
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
    size_t j;
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

float synQuantSimple::debiasedFgBgCompare(unsigned debiasMethod){

}