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
synQuantSimple::synQuantSimple(Mat *_srcVolume, float _src_var, segParameter &p4segVol, odStatsParameter &p4odStats){
    srcVolumeUint8 = _srcVolume;
    imArray = (unsigned char*)srcVolumeUint8->data;
    src_var = _src_var;

    componentTree3d(p4segVol, p4odStats);
    float * zscore_pointer = & (*zscore_list.begin()); // a pointer to the address of zscore_list.begin()
    zMap = new Mat(srcVolumeUint8->dims, srcVolumeUint8->size, CV_32F, zscore_pointer); //*zscore_list.begin()
    cell_num = floatMap2idMap(zMap, *idMap, 26);
}
synQuantSimple::synQuantSimple(singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats){
    srcVolumeUint8 = &seed.volUint8;
    imArray = (unsigned char*)srcVolumeUint8->data;
    vector<float> var_vals_in_fg = extractValsGivenMask(&seed.varMap, &seed.fgMap, CV_8U, 0);
    src_var = vec_mean(var_vals_in_fg);
    // 1. find best threshold and apply some basic morphlogical operation (ordstat4fg.m)
    // PS: the region from this step will has no overlap with other cell's terriotory
    // PPS: the cell territory may contain more than one cell, we need verify it using gaps.
    cellTerritoryExtractFromSeed(seed, p4odStats);
    refineCellTerritoryWithSeedRegion(seed, p4segVol);

    // 2. update seed's score map based on idMap (fg)
    normalize(seed.eigMap2d, seed.score2d, 0.001, 1, NORM_MINMAX, CV_32F, *idMap);
    normalize(seed.eigMap3d, seed.score3d, 0.001, 1, NORM_MINMAX, CV_32F, *idMap);

    // 3. segment idMap (fg) based on gaps from principal curvature
    gapBasedRegionSegment(seed, p4segVol, p4odStats);
    // 4. refine the output from step 3 based on size and other prior knowledge

    // 5. test if the idMap (fg) is too small, if so, enlarge it and re-do the previous steps

}
/**
 * @brief refineCellsTerritoryWithSeedRegions: the last time of refinement based size and shrinking
 * @param seed
 * @param p4segVol
 */
void synQuantSimple::refineCellsTerritoriesWithSeedRegions(singleCellSeed &seed, segParameter &p4segVol){

}
/**
 * @brief refineFgWithSeedRegion, extract the largest region if fg and remove those
 * voxels that has no z-direction neighbors
 * @param seed
 * @param p4segVol
 */
void synQuantSimple::refineCellTerritoryWithSeedRegion(singleCellSeed &seed, segParameter &p4segVol){
    Mat labelMap;
    int n = connectedComponents3d(idMap, labelMap, 6);
    if (n > 1){
        int id = largestRegionIdExtract(&labelMap, n, &seed.seedMap);
        *idMap = labelMap == id;
    }
    // remove pixels that happen only at one slice
    Mat tmp_map;
    idMap->copyTo(tmp_map);
    FOREACH_ijk_ptrMAT(idMap){
        if(tmp_map.at<unsigned char>(i,j,k) > 0){
            if((k-1>=0 && tmp_map.at<unsigned char>(i,j,k-1) == 0) &&
                (k+1<tmp_map.size[2] && tmp_map.at<unsigned char>(i,j,k+1) == 0)){
                idMap->at<unsigned char>(i,j,k) = 0;
            }
        }
    }
    int radius[] = {1,1,0};
    volumeDilate(idMap, tmp_map, radius, MORPH_ELLIPSE);
    n = connectedComponents3d(&tmp_map, labelMap, 26);
    removeSmallCC(labelMap, n, p4segVol.min_cell_sz, false);
    *idMap = labelMap > 0;
}

/**
 * @brief fgGapRemoval: remove the gap defined by principal curvature. We jointly consider the 2d and 3d
 * principal curvature, since 3d may be too liberal to remove too much areas
 * @param seed
 * @param p4segVol
 */
void synQuantSimple::fgGapRemoval(singleCellSeed &seed, segParameter &p4segVol){
    bitwise_and(*idMap, seed.gap3dMap == 0, seed.fgMapGapRemoved);

    if (!isempty(&seed.fgMapGapRemoved, CV_8U)){
        Mat seedMapFrom2dMap, mapUnion, mapUnion_label;
        bitwise_and(*idMap, seed.gap2dMap == 0, seedMapFrom2dMap);// or use -
        bitwise_or(seedMapFrom2dMap, seed.fgMapGapRemoved, mapUnion); // or use +
        int numCC = connectedComponents3d(&mapUnion, mapUnion_label, p4segVol.connect4fgGapRemoval);
        Mat newSeedMap; // CV_8U
        bool found = findUnrelatedCC(&mapUnion_label, numCC, &seed.fgMapGapRemoved, newSeedMap);
        if(found){
            seed.fgMapGapRemoved = seed.fgMapGapRemoved + newSeedMap; //CV_8U
        }
    }
}
/**
 * @brief gapBasedRegionSegment: use principal curvature to test if current fg contains > 1 cells
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void synQuantSimple::gapBasedRegionSegment(singleCellSeed &seed, segParameter &p4segVol, odStatsParameter &p4odStats){
    // 1. remove gaps defined by principal curvature
    fgGapRemoval(seed, p4segVol);
    // 2. check if gaps divided the region into multiple valid seeds
    Mat label_map;
    int n = connectedComponents3d(&seed.fgMapGapRemoved, label_map, p4segVol.neiMap);
    removeSmallCC(label_map, n, p4segVol.min_seed_size, REARRANGE_IDS);
    if(n <= 1){
        return;
    }
    // 3. gap test: test if the gap is true, if so, split the fg with respect to the multiple seeds
    gapTest2SplitCellTerritory(&label_map, n, seed, p4segVol, p4odStats);
}
/**
 * @brief gapTest2SplitCellTerritory: Given a seeds_map, test if the gap between any two seeds are true or false.
 * If true, keep them. Otherwise, merge the two seeds divided by the gap. If the gap is indeed dimmer than the
 * two seeds it split, it is true. Otherwise false.
 * @param seeds_Map:cv_32s
 * @param n
 * @param seed
 * @param p4segVol
 * @param p4odStats
 */
void synQuantSimple::gapTest2SplitCellTerritory(Mat* seeds_Map /*CV_32S*/, int n,
                                                singleCellSeed &seed, segParameter &p4segVol,
                                                odStatsParameter &p4odStats){
    // 1. first grow all seeds until they touch with each other
    Mat grown_seedMap;
    Mat scoreMap = seed.score2d + seed.score3d;
    bool link_bg2sink = false; // this can force seeds to grow as much as they can
    regionGrow(seeds_Map, n, grown_seedMap, &scoreMap, &seed.fgMap,
               p4segVol.growConnectInRefine, p4segVol.graph_cost_design, link_bg2sink);

    // 2. for the connected seeds, test if the gap among them are true
    // the gap are choosen by dilating the seed regions
    vector<bool> gap_tested_true(n*n, false);
    float p_treshold = p4odStats.gapTestThreshold / (p4segVol.gapTestMinMaxRadius[1] - p4segVol.gapTestMinMaxRadius[0]);
    vector<size_t> real_gap_idx;
    for(int r = p4segVol.gapTestMinMaxRadius[0]; r <= p4segVol.gapTestMinMaxRadius[1]; r++){
        vector<vector<size_t>> gap_idx_list;
        extractGapVoxel(&grown_seedMap, &seed.fgMap, n, r, gap_idx_list, gap_tested_true);
        FOREACH_i(gap_idx_list){
            if(gap_idx_list[i].size() > 0){
                int id0 = i / n + 1;
                int id1 = i % n + 1;
                vector<vector<size_t>> nei_list0, nei_list1;
                Mat cur_fg = grown_seedMap == id0;
                neighbor_idx_2d(gap_idx_list[i], &cur_fg, nei_list0, 2*r+1);
                cur_fg = grown_seedMap == id1;
                neighbor_idx_2d(gap_idx_list[i], &cur_fg, nei_list1, 2*r+1);

                double var_sum = 0;
                size_t var_idx_num = 0;
                vector<float> r0, r1, skipped, gap;

                gap = extractValsGivenIdx(&seed.volStblizedFloat, gap_idx_list[i], CV_32F);
                var_sum += extractSumGivenIdx(&seed.stblizedVarMap, gap_idx_list[i], CV_32F);
                var_idx_num += gap_idx_list[i].size();
                for(size_t j = p4odStats.gapTestSkippedBandWidth; j < nei_list0.size(); j++){
                    vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list0[j], CV_32F);
                    r0.insert(r0.end(), curr_vals.begin(), curr_vals.end());
                    var_sum += extractSumGivenIdx(&seed.stblizedVarMap, nei_list0[j], CV_32F);
                    var_idx_num += nei_list0[j].size();
                }
                for(size_t j = p4odStats.gapTestSkippedBandWidth; j < nei_list1.size(); j++){
                    vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list1[j], CV_32F);
                    r1.insert(r1.end(), curr_vals.begin(), curr_vals.end());
                    var_sum += extractSumGivenIdx(&seed.stblizedVarMap, nei_list1[j], CV_32F);
                    var_idx_num += nei_list1[j].size();
                }
                float cur_std = (float)sqrt(var_sum / var_idx_num);
                float mu, sigma;
                float p0, p1;
                if(p4odStats.gapTestMethod == GAP_LOCALORDERSTATS){
                    float sum_stats0 = vec_mean(r0) - vec_mean(gap);
                    float sum_stats1 = vec_mean(r1) - vec_mean(gap);
                    gap.resize(0);
                    for(size_t j = 0; j < p4odStats.gapTestSkippedBandWidth; j++){
                        vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list0[j], CV_32F);
                        skipped.insert(skipped.end(), curr_vals.begin(), curr_vals.end());
                    }
                    for(size_t j = 0; j < p4odStats.gapTestSkippedBandWidth; j++){
                        vector<float> curr_vals = extractValsGivenIdx(&seed.volStblizedFloat, nei_list1[j], CV_32F);
                        skipped.insert(skipped.end(), curr_vals.begin(), curr_vals.end());
                    }
                    orderStatsKSection(r0, gap, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p0 = zscore2palue((sum_stats0 - mu) / sigma);
                    p1 = zscore2palue((sum_stats1 - mu) / sigma);
                }else if(p4odStats.gapTestMethod == GAP_ORDERSTATS){
                    float sum_stats0 = vec_mean(r0) - vec_mean(gap);
                    float sum_stats1 = vec_mean(r1) - vec_mean(gap);

                    orderStatsKSection(r0, gap, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p0 = zscore2palue((sum_stats0 - mu) / sigma);

                    orderStatsKSection(r1, gap, skipped, mu, sigma);
                    mu *= cur_std;
                    sigma *= cur_std;
                    p1 = zscore2palue((sum_stats1 - mu) / sigma);
                }

                if (p0 <= p_treshold && p1 <= p_treshold){ // gap is true
                    gap_tested_true[i] = true;
                    real_gap_idx.insert(real_gap_idx.end(), gap_idx_list[i].begin(), gap_idx_list[i].end());
                }
            }
        }
        gap_idx_list.clear();
    }
    vector<vector<int>> groups;
    FOREACH_i(gap_tested_true){
        if(gap_tested_true[i]){
            int target0 = i / n + 1;
            int target1 = i % n + 1;
            bool found_group = false;
            for(size_t j = 0; j < groups.size(); j++){
                for(size_t k = 0; k < groups[j].size(); k++){
                    if(target0 == groups[j][k]){
                        groups[j].push_back(target1);
                        found_group = true;
                        break;
                    }else if(target1 == groups[j][k]){
                        groups[j].push_back(target0);
                        found_group = true;
                        break;
                    }
                }
                if(found_group) break;
            }
            if(!found_group){
                groups.push_back({target0, target1});
            }
        }
    }
    vector<int> label_map (n);
    FOREACH_i(label_map){
        label_map[i] = i + 1;
    }
    FOREACH_i(groups){
        //vec_unique(groups[i]); //sorted also
        for(size_t j = 1; j < groups[i].size(); j++){
            assert(label_map[groups[i][j] - 1] == groups[i][j] || label_map[groups[i][j] - 1] == groups[i][0]);// "Duplicated Regions in two groups");
            label_map[groups[i][j] - 1] = groups[i][0];
        }
    }
    FOREACH_i_MAT(grown_seedMap){
        if(grown_seedMap.at<int>(i) > 0){
            grown_seedMap.at<int>(i) = label_map[grown_seedMap.at<int>(i) - 1];
        }
    }
    vector<size_t> non_used;
    rearrangeIdMap(&grown_seedMap, *idMap, non_used);
    FOREACH_i(real_gap_idx){ // remove the real gaps
        idMap->at<int>(real_gap_idx[i]) = 0;
    }
}

/**
 * @brief synQuantSimple::cellExtractFromSeed, currently, we use the exact version. But this function should also be able
 * to be implemented using component tree condition on the assumption that the intensity has such trend:
 * inner region > boundary > neighbor. The boundary selection can be fulfilled with the similar way as neighbor selection.
 * @param seed
 * @param p4odStats
 */
void synQuantSimple::cellTerritoryExtractFromSeed(singleCellSeed &seed, odStatsParameter &p4odStats){
    //Scalar seed_val = mean(seed.volUint8, seed.fgMap);// scalar is a vector<double> with length 4 for RGBA data
    int ub = MAX(30, round(mat_mean(&seed.volUint8, CV_8U, seed.idx_yxz_cropped)));
    int lb = MAX(5, round(mean(seed.volUint8, seed.fgMap)[0]));
    //Mat curr_seed_map = seed.idMap == seed.id; //replaced by seed.seedMap
    if (ub <= lb){
        *idMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
        return;
    }
    vector<float> zscore (ub-lb+1);

    Mat otherCellTerritory, valid_nei, cur_reg, cur_valid_nei; // uint8
    int shift_othercell[] = {3,3,1};
    volumeDilate(&seed.otherIdMap, otherCellTerritory, shift_othercell, MORPH_ELLIPSE);
    bitwise_and(otherCellTerritory, seed.seedMap==0, otherCellTerritory); // remove the area covered by seed region
    bitwise_and(otherCellTerritory == 0, seed.fgMap, valid_nei); // remove the area not inside fgMap
    double seed_sz = (double)seed.idx_yxz.size();
    max_exist_zscore = -1;
    maxZ_intensity_level = -1;
    FOREACH_i(zscore){
        bitwise_and(seed.volUint8 >= (ub-i), otherCellTerritory == 0, cur_reg);
        validSingleRegionExtract(cur_reg, &seed.seedMap, p4odStats.connectInSeedRefine);
        bitwise_and(valid_nei, cur_reg==0, cur_valid_nei);

        size_t reg_sz = overlap_mat_vec(&cur_reg, CV_8U, seed.idx_yxz_cropped, 0);
        if ((reg_sz / seed_sz) < 0.5){
            continue;
        }
        float cur_zscore = debiasedFgBgBandCompare(&cur_reg, &valid_nei, &seed, p4odStats);
        if (max_exist_zscore < cur_zscore){
            maxZ_intensity_level = (ub-i);
            max_exist_zscore = cur_zscore;
        }
    }
    if (max_exist_zscore > 0){
        bitwise_and(seed.volUint8 >= maxZ_intensity_level, otherCellTerritory == 0, *idMap);
        validSingleRegionExtract(*idMap, &seed.seedMap, p4odStats.connectInSeedRefine);
    }else{
        *idMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
        return;
    }
    // refine the resultant region
    if (isempty_mat_vec(idMap, CV_8U, seed.idx_yxz_cropped, 0)){
        *idMap = Mat::zeros(seed.volUint8.dims, seed.volUint8.size, CV_8U);
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
    int dilate_size[] = {p4odStats.gap4fgbgCompare, p4odStats.gap4fgbgCompare, 0};
    volumeDilate(cur_reg, cur_reg_dilate, dilate_size, MORPH_ELLIPSE);
    //volumeErode(cur_reg, cur_reg_erode, dilate_size, MORPH_ELLIPSE);

    Mat fg_center, fg_band, bg_band;
    //bitwise_and(*cur_reg, cur_reg_erode == 0, fg_band);
    bitwise_and(*validNei, cur_reg_dilate, bg_band);
    volumeDilate(&bg_band, fg_band, dilate_size, MORPH_ELLIPSE);
    bitwise_and(fg_band, *cur_reg, fg_band);
    fg_center = seed->fgMap - fg_band - bg_band; // key here
    vector<float> fg_band_vals, bg_band_vals, fg_center_vals;
    fg_band_vals = extractValsGivenMask(&seed->volUint8, &fg_band, CV_8U, 0);
    bg_band_vals = extractValsGivenMask(&seed->volUint8, &bg_band, CV_8U, 0);
    fg_center_vals = extractValsGivenMask(&seed->volUint8, &fg_center, CV_8U, 0);

    float mu, sigma, zscore = 0;
    if (p4odStats.fgSignificanceTestWay == KSEC){
        orderStatsKSection(fg_band_vals, bg_band_vals, fg_center_vals, mu, sigma);
        float sum_stats = vec_mean(fg_band_vals) - vec_mean(bg_band_vals);
        zscore = (sum_stats - mu*sqrt(src_var))/(sigma*sqrt(src_var));
    }else{
        //TODO: other ways haven't been implemented here
    }
    return zscore;
}
