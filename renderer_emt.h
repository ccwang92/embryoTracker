#ifndef RENDERER_EMT_H
#define RENDERER_EMT_H


#include <GL/glu.h>
//#include "src_3rd/3drenderer/GLee2glew.h"
//#include <glew/GL/glew.h>
#include "src_3rd/basic_c_fun/color_xyz.h"
#include "data_importer.h"
//#include "src_3rd/3drenderer/glsl_r.h"
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL related
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//@change from inline to define: 2020-2-10 RZC

#define supported_MultiSample() \
( \
    GLeeInit() || (GLEE_ARB_multisample>0) \
)
#define supported_Tex3D() \
( \
    GLeeInit() || (GLEE_EXT_texture3D>0 || GLEE_VERSION_2_0>0) \
)
//#define supported_Tex3D()  (getMaxTexSize3D()>0) //// (VERSION_2_0==0 when ssh -X, EXT_texture3D==0 when NVIDIA)
#define supported_TexCompression() \
( \
    GLeeInit() || (GLEE_ARB_texture_compression>0) \
)
#define supported_TexNPT() \
( \
    GLeeInit() || (GLEE_ARB_texture_non_power_of_two>0) \
)
#define supported_GL2() \
( \
    GLeeInit() || (GLEE_VERSION_2_0>0) \
)
#define supported_GLSL() \
( \
    GLeeInit() || (GLEE_ARB_shading_language_100>0 || GLEE_VERSION_2_0>0) \
)
#define supported_PBO() \
( \
    GLeeInit() || (GLEE_ARB_pixel_buffer_object>0) \
)

/////////////////////////////////////////////////////////////////////////////

#define CHECK_GLErrorId_throw() \
{ \
    GLenum id = glGetError(); \
    switch (id) \
    { \
    case GL_NO_ERROR: \
        break; \
    case GL_INVALID_ENUM: \
    case GL_INVALID_VALUE: \
    case GL_INVALID_OPERATION: \
        break; \
    case GL_STACK_OVERFLOW: \
    case GL_STACK_UNDERFLOW: \
    case GL_OUT_OF_MEMORY: \
    case GL_TABLE_TOO_LARGE: \
        throw id; break; \
    default: \
        throw id; \
    } \
}
#define CHECK_GLErrorString_throw() \
try { \
    CHECK_GLErrorId_throw(); \
} catch(GLenum id) { \
    const char* str = (const char*)gluErrorString(id); \
    throw str; \
}
#define CHECK_GLError_print() CheckGLError_print(__FILE__, __LINE__)
#define Q_CSTR(qs)  ( (qs).toStdString().c_str() )
inline int CheckGLError_print(const char *file, int line)
{
    GLenum glErr;
    while ((glErr = glGetError()) != GL_NO_ERROR)
    {
        const GLubyte* sError = gluErrorString(glErr);
        if (sError)
        std::cerr << "GL Error #" << glErr << "(" << sError << ") " << " in file(" << file << ") at line(" << line << ")\n";
        else
        std::cerr << "GL Error #" << glErr << " (no message available)" << " in file(" << file << ") at line(" << line << ")\n";
    }
    return glErr;
}
#define MAX_3DTEX_SIZE  256
// if error then close
// clean memory before MessageBox, otherwise MessageBox maybe could not be created correctly
#define ERROR_MessageBox(title, type, what) { \
    b_error = true; \
    cleanData(); \
    QMessageBox::critical( 0, title, QObject::tr("%1: OUT OF MEMORY or GL ERROR.\n---%2 exception: %3")\
            .arg(title).arg(type).arg(what) + "\n\n" + \
        QObject::tr("3D View: Please close some images or views to release memory, then try again.\n\n") ); \
}
#define CATCH_handler( func_name ) \
    catch (std::exception & e) { \
        \
        qDebug("   catch: data4dp  = @%0p", data4dp); \
        qDebug("   catch: rgbaBuf = @%0p", rgbaBuf); \
        \
        ERROR_MessageBox( func_name, "std", e.what() ); \
        return; \
        \
    } catch (const char* s) { \
        \
        ERROR_MessageBox( func_name, "GL", s ); \
        return; \
        \
    } catch (...) { \
        \
        ERROR_MessageBox( func_name, "UNKOWN", "unknown exception" ); \
        return; \
        \
    }

#define MESSAGE_ASSERT(s) \
{\
    if (!(s)) \
        QMessageBox::critical(0, "ASSERT", QObject::tr("ASSERT(%1) in file(%2) at line(%3)").arg(#s).arg(__FILE__).arg(__LINE__)); \
    Q_ASSERT(s); \
}

#define DELETE_AND_ZERO(p)	{ if ((p)!=NULL) delete (p); (p) = NULL; }
//special definition: there are some cases these following definitions are missing, and always define here if not defined before. 090730
// for mac 10.4 __DARWIN__ ?
#ifndef GL_SHADING_LANGUAGE_VERSION_ARB
#define GL_SHADING_LANGUAGE_VERSION_ARB  0x8B8C
#endif
#ifndef GL_MAX_DRAW_BUFFERS_ARB
#define GL_MAX_DRAW_BUFFERS_ARB  0x8824
#endif
#ifndef GL_MAX_COLOR_ATTACHMENTS_EXT
#define GL_MAX_COLOR_ATTACHMENTS_EXT  0x8CDF
#endif

#define N_CHANNEL 3
#define FILL_CHANNEL 4  // power_of_two

#define IF_VOL_SHADER  (tryVolShader && shader && !b_selecting)
#define IF_OBJ_SHADER  (tryObjShader && shader && !b_selecting)

enum v3dr_DataClass { dcDataNone=0,
                dcVolume=1,
                dcSurface=2
                };
enum v3dr_VolSlice { vsSliceNone=0,
                vsXslice=1,
                vsYslice=2,
                vsZslice=3,
                vsFslice=4,
                };
enum v3dr_SurfaceType { stSurfaceNone=0,
                stImageMarker=1,
                stLabelSurface=2,
                stNeuronSegment=3,
                stNeuronStructure=4,
                stPointCloud=5,
                stPointSet=6,
                };
enum RenderMode {rmCrossSection=0,
                rmAlphaBlendingProjection,
                rmMaxIntensityProjection,
                rmMinIntensityProjection,
                };
class Renderer_EmT
{
    //friend class V3dr_colormapDialog; // for access colormap, colormap_curve

public:
    Renderer_EmT(void* widget);
    virtual ~Renderer_EmT();
    //virtual const int class_version() {return 2;}
    void makeCurrent();
    void blendBrighten(float fbright, float fcontrast);
// link to Data
    virtual void setupData(void* data);
    virtual void cleanData();                      // makeCurrent
    virtual const bool has_image()   {return (size4>0);}
    virtual const BoundingBox getDataBox() {return dataBox;}
    virtual void getLimitedDataSize(int size[5]) {for (int i=0;i<5;i++) size[i]=bufSize[i];};
// link to OpenGL window
    virtual void initialize();
    virtual void setupView(int width, int height);		//link to QGLWidget::resizeGL
    virtual void setProjection(); 			// called by setupView & selectObj
    virtual void paint();
// link to View control
    virtual void toggleShader();
    virtual void toggleObjShader();
    virtual void applyColormapToImage();
    virtual RGB8 lookupColormap(RGB8 inC, int op);
    virtual void toggleTexStream();

//	virtual void togglePolygonMode() {polygonMode = (polygonMode +1) %5;} // FILL,LINE,POINT, transparent,outline

// link to Rendering function
//protected:
    virtual void loadObj()	{boundingBox=surfBoundingBox=UNIT_BoundingBox;}	// MUST makeCurrent for concurrent contexts, 081025
    virtual void loadVol();
    virtual void cleanVol();

    virtual void loadShader();  // called by initialize()  	// makeCurrent
    virtual void cleanShader(); // called by ~Renderer_EmT 	// makeCurrent

    virtual void prepareVol();
    virtual void renderVol();
    virtual void drawVol();
    virtual void drawObj();  // called by paint()
    //virtual void drawVol();  // called by paint() //use default.

    virtual void equAlphaBlendingProjection();
    virtual void equMaxIntensityProjection();
    virtual void equMinIntensityProjection();
    virtual void equCrossSection();
    virtual void shaderTexBegin(bool stream); //090722 because different shader for XYZ/F-slice
    virtual void shaderTexEnd();

    virtual void initColormap();

    virtual bool supported_TexStream(); // called by loadVol()
    virtual void setupTexStreamBuffer(); // called by subloadTex()
    virtual void cleanTexStreamBuffer(); // called by ~Renderer_EmT 	// makeCurrent
    virtual void setupStackTexture(bool bfirst);
    virtual void _streamTex(int stack_i, int slice_i, int step, int slice0, int slice1);
    virtual void _streamTex_end();
    virtual bool _streamTex_ready();
// rendering execution function from renderer_gl1
    virtual void setObjLighting();
    virtual void disObjLighting();
    virtual void beginHighlight();
    virtual void endHighlight();

    virtual void setUnitVolumeSpace();
    virtual void drawUnitVolume();
//rendering execution function from renderer
    virtual void setViewClip(float f);
    virtual void enableViewClipPlane();
    virtual void disableViewClipPlane();
//    virtual void setupStackTexture(bool bfirst);
    virtual int  _getBufFillSize(int w);
    virtual int  _getTexFillSize(int w);
    virtual void setTexParam3D();
    virtual void setTexParam2D();
    virtual void drawStackZ(float direction, int section, bool t3d, bool stream);
    virtual void drawStackY(float direction, int section, bool t3d, bool stream);
    virtual void drawStackX(float direction, int section, bool t3d, bool stream);
    virtual void _drawStack( double ts, double th, double tw,
            double s0, double s1, double h0, double h1, double w0, double w1,
            double ds, int slice0, int slice1, int thickness,
            GLuint tex3D, GLuint texs[], int stack_i,
            float direction, int section, bool t3d, bool stream);

//    virtual bool supported_TexStream()								{return false;} //091003
//    virtual void setupTexStreamBuffer()							{tex_stream_buffer = false;}; //091003
//    virtual void cleanTexStreamBuffer()							{tex_stream_buffer = false;}; //091003
//    virtual void _streamTex(int stack_i, int slice_i, int step, int slice0, int slice1) {}; //091003
//    virtual void _streamTex_end()													 	{}; //091003
//    virtual bool _streamTex_ready()									{return false;} //091015

    virtual void setupFrontSliceBuffer();
    virtual void drawUnitFrontSlice(int line=0);

    virtual void drawBackFillVolCube();
    virtual void drawCrossLine(float lineWidth=1);
// link to rendering functions from renderer_gl1
    virtual void subloadTex(V3DLONG timepoint, bool bfisrt=false);	// called by loadVol, drawVol
/*	virtual void equAlphaBlendingProjection();			// for blending equation
    virtual void equMaxIntensityProjection();	// for blending equation
    virtual void equMinIntensityProjection();	// for blending equation
    virtual void equCrossSection();				// for blending equation
    virtual void shaderTexBegin(bool stream) 	{};		// for texture2D/3D shader
    virtual void shaderTexEnd()	 				{};	*/	// for texture2D/3D shader
    virtual void updateVolCutRange();
    virtual void updateBoundingBox();
    virtual void updateThicknessBox();

public: // renderer
    RenderMode renderMode;
    int sShowTrack, curChannel;
     int sShowRubberBand; // ZJL 1109221

    bool bShowBoundingBox, bShowBoundingBox2, bShowAxes, bOrthoView;
    // TDP 201601 - provide optional XY translation arrows directly in 3D view to navigate to adjacent/overlapping ROI
    bool bShowXYTranslateArrows;
    int iPosXTranslateArrowEnabled, iNegXTranslateArrowEnabled, iPosYTranslateArrowEnabled, iNegYTranslateArrowEnabled;
    BoundingBox* posXTranslateBB;
    BoundingBox* negXTranslateBB;
    BoundingBox* posYTranslateBB;
    BoundingBox* negYTranslateBB;
    bool bShowCSline, bShowFSline, bFSlice, bXSlice, bYSlice, bZSlice;
    float CSbeta, alpha_threshold;
    RGBA32f color_background, color_background2, color_line, color_proxy;

    int sShowMarkers, sShowSurfObjects, markerSize;
    bool b_showMarkerLabel, b_showMarkerName, b_showCellName, b_surfStretch, b_surfZLock;

    int lineType, lineWidth, nodeSize, rootSize;
    int polygonMode, tryObjShader;
    int tryTexNPT, tryTex3D, tryTexCompress, tryVolShader, tryTexStream;
    const char* try_vol_state();
    int editinput;
    void rgba3d_r2gray(RGBA8* rgbaBuf, V3DLONG bufSize[5]);
    void data4dp_to_rgba3d(Image4DProxy<Image4DSimple>& img4dp, V3DLONG dim5,
                           V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
                           V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
                           RGBA8* rgbaBuf, V3DLONG bufSize[5]);
    void data4dp_to_rgba3d(unsigned char* data4dp, V3DLONG dim1, V3DLONG dim2, V3DLONG dim3, V3DLONG dim4, V3DLONG dim5,
                           V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
                           V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
                           RGBA8* rgbaBuf, V3DLONG bufSize[5]);
// internal state
//protected:
    int volTimePoint, volTimeOffset;
    BoundingBox boundingBox, surfBoundingBox;
    double thickness; //changed from int to double, PHC, 090215
    V3DLONG xCut0,xCut1, yCut0,yCut1, zCut0,zCut1;            // for volume
    float xClip0,xClip1, yClip0,yClip1, zClip0,zClip1;    // for surface
    float viewClip;

    V3DLONG screenW, screenH;
    float viewDistance, viewNear, viewFar, viewAngle, zoomRatio;

    bool b_error;
    bool b_selecting;

    bool b_limitedsize;
    float *depth_buffer;
    RGBA8 *total_rgbaBuf, *rgbaBuf;
    float sampleScale[5];
    V3DLONG bufSize[5]; //(x,y,z,c,t) 090731: add time dim

    XYZ curveStartMarker; // ZJL
    int neuronColorMode;
    int dispConfLevel;

public: // renderer_gl1
    void* widget; // get the gl widget
    void* _idep;

    int data_unitbytes;
    unsigned char* data4dp;
    unsigned char**** data4d_uint8;
    // data4d_uint8[dim4][dim3][dim2][dim1]
    V3DLONG dim1, dim2, dim3, dim4, dim5;
    V3DLONG start1, start2, start3, start4, start5;
    V3DLONG size1, size2, size3, size4, size5;
    BoundingBox dataBox;
    BoundingBox dataViewProcBox; //current clip box that data are visible (and thus are processable). 091113 PHC

    bool texture_unit0_3D, tex_stream_buffer, drawing_fslice;
    GLenum texture_format, image_format, image_type;
    GLuint tex3D, texFslice;
    GLuint *Ztex_list, *Ytex_list, *Xtex_list;
    RGBA8 *Zslice_data, *Yslice_data, *Xslice_data, *Fslice_data;
    RGBA8 *rgbaBuf_Yzx, *rgbaBuf_Xzy;
    float thicknessX, thicknessY, thicknessZ;
    float sampleScaleX, sampleScaleY, sampleScaleZ;
    int imageX, imageY, imageZ, imageT;
    int safeX, safeY, safeZ;
    int realX, realY, realZ, realF;
    int fillX, fillY, fillZ, fillF;
    GLdouble volumeViewMatrix[16]; // for choosing stack direction
    float VOL_X1, VOL_X0, VOL_Y1, VOL_Y0, VOL_Z1, VOL_Z0;
    int VOLUME_FILTER;
    RGBA32f SLICE_COLOR; // proxy geometry color+alpha
     bool b_renderTextureLast;
    double currentTraceType;
    bool useCurrentTraceTypeForRetyping;

    RGBA8 currentMarkerColor;//added by ZZ 05142018

    float zThick;

public: // renderer_gl2
    RGBA8 colormap[FILL_CHANNEL][256];      // [n-channel][256-intensity]
    cwc::glShaderManager SMgr;
    cwc::glShader *shader, *shaderTex2D, *shaderTex3D, *shaderObj;

    GLuint texColormap; // nearest filter, [x-coord for intensity][y-coord for channel]
    // RGBA8 colormap[FILL_CHANNEL][256];      // [n-channel][256-intensity]
    QPolygonF colormap_curve[N_CHANNEL][4]; // [n-channel][RGBA]

    GLuint pboZ, pboY, pboX;
    GLenum pbo_texture_format, pbo_image_format, pbo_image_type;

private:
    void init_members()
    {
        /********Renderer*******/
        //b_useClipBoxforSubjectObjs = true;

        sShowTrack = 0;
        curChannel = 0;
          sShowRubberBand = 0; // ZJL 110921

          curveStartMarker = XYZ(0,0,0);

        bShowBoundingBox = true;
        bShowBoundingBox2 = false;
        bShowAxes = true;
        bOrthoView = false;

        bShowXYTranslateArrows = 0;
        iPosXTranslateArrowEnabled = 0;
        iNegXTranslateArrowEnabled = 0;
        iPosYTranslateArrowEnabled = 0;
        iNegYTranslateArrowEnabled = 0;
        posXTranslateBB = 0;
        negXTranslateBB = 0;
        posYTranslateBB = 0;
        negYTranslateBB = 0;

        bShowCSline = true;
        bShowFSline = true;
        bXSlice = bYSlice = bZSlice = true;
        bFSlice = true; //100729 RZC: it's not slow, default turn on.
        CSbeta = alpha_threshold = 0;
        thickness = 1;

        color_background = XYZW(.1f, .1f, .25f, 1); // background color for volume
        color_background2 = XYZW(.8f, .85f, .9f, 1); // background color for geometric object only
        color_line = XYZW(.5f,.5f,.7f, 1);
        color_proxy = XYZW(1, 1, 1, 1);

        sShowMarkers = 2;
        sShowSurfObjects = 2;
        markerSize = 10;
        b_showMarkerLabel = true;
        b_showMarkerName = false; // by Lei Qu, 110425
        b_surfStretch = true;
        b_showCellName    = false;

        lineType = 1; lineWidth = 1; nodeSize = 0, rootSize = 5;
        polygonMode = 0;
        tryObjShader = 0;

        tryTexNPT = 0;
        tryTex3D = 0;
        tryTexCompress = 1;
        tryTexStream = 1;
        tryVolShader = 1;

        volTimePoint = volTimeOffset = 0;

        //protected -------------------------------------

        xCut0 = yCut0 = zCut0 = -1000000; // no cut
        xCut1 = yCut1 = zCut1 = 1000000;  // no cut
        xClip0 = yClip0 = zClip0 = -1000000; // no clip
        xClip1 = yClip1 = zClip1 = 1000000;  // no clip
        viewClip = 1000000;  // no clip

        //// perspective view frustum
        screenW = screenH = 0;
        viewAngle = 31;
        zoomRatio = 1;
        viewNear = 1;
        viewFar = 10;
        viewDistance = 5;

        b_error = b_selecting = false;

        b_limitedsize = false;
        depth_buffer = 0;
        total_rgbaBuf = rgbaBuf = 0;
        for (int i=0; i<5; i++)
        {
            sampleScale[i]=1; bufSize[i]=0;
        }


        renderMode = rmMaxIntensityProjection;
        //selectMode = smObject;

        //refineMode = smCurveRefine_fm;

        //ui3dviewMode = Vaa3d;
        editinput=0;
        neuronColorMode=0;
        dispConfLevel=INT_MAX;
        /*******Renderer gl1***********/
//        lastSliceType = vsSliceNone;
//        currentMarkerName = -1;
//        curEditingNeuron = -1;
//        realCurEditingNeuron_inNeuronTree = -1;

//        highlightedNode = -1; //Added by ZMS 20151203 highlight initial node we are going to extend.
//        highlightedEndNode = -1; //Added by ZMS 20151203 highlight final node we are going to extend.
//        selectedStartNode = -1;
//        highlightedEndNodeChanged = false;
//        rotateAxisBeginNode = XYZ(0, 0, 1);
//        rotateAxisEndNode = XYZ(0, 0, 0);

//        childHighlightMode = false;
//        showingGrid = false;
        _idep=0;
//        isSimulatedData=false;
        data_unitbytes=0;
        data4dp = 0;
        data4d_uint8 = 0;
        dim1=dim2=dim3=dim4=dim5=0;
        start1=start2=start3=start4=start5=0;
        size1=size2=size3=size4=size5=0;

        texture_unit0_3D = tex_stream_buffer = drawing_fslice = false;
        texture_format = image_format = image_type = -1;
        tex3D = texFslice= 0;
        Ztex_list = Ytex_list = Xtex_list = 0;
        Zslice_data = Yslice_data = Xslice_data = Fslice_data = 0;
        rgbaBuf_Yzx = rgbaBuf_Xzy = 0;
        thicknessX = thicknessY = 1; thicknessZ = 1;
        sampleScaleX = sampleScaleY = sampleScaleZ = 1;
        imageX = imageY = imageZ = imageT = 0;
        safeX = safeY = safeZ = 0;
        fillX = fillY = fillZ = fillF = 0;
        realX = realY = realZ = realF = 0;
        VOL_X1 = VOL_Y1 = VOL_Z1 = 1;
        VOL_X0 = VOL_Y0 = VOL_Z0 = 0;
        VOLUME_FILTER = 1;
        SLICE_COLOR = XYZW(1,1,1,1);

//        b_grabhighrez = false; //120717, PHC
//        b_imaging = false; //101008
//        b_ablation = false; //120506
//        b_lineAblation = false;
          b_renderTextureLast = false;
          /*edit_seg_id = -1; // ZJL 110913
          draggedCenterIndex = -1; // ZJL 110921
          nDragWinSize = 3; // better be odd, ZJL 110921
          bInitDragPoints = false;
          // for curve testing
          bTestCurveBegin=false;

          b_editDroppedNeuron = false*/; //20150527, PHC

//          highlightedNodeType = -1; //20170804 RZC
          currentTraceType=3;
          useCurrentTraceTypeForRetyping = false;
//        cuttingZ = false;
//        cuttingXYZ = false;
//        zMin =-1.0;
//        zMax = 1.0;
//        initColorMaps();
//        gridSpacing = 10.0;
        currentMarkerColor.r=255;
        currentMarkerColor.g=0;
        currentMarkerColor.b=0;
//        deleteKey = 0;
        /******Renderer gl2******/
        shader = shaderTex2D=shaderTex3D = shaderObj = 0;
        texColormap = 0;
        pboZ = pboY = pboX = 0;
        pbo_texture_format = pbo_image_format = pbo_image_type = -1;
    }

// helping function

inline void set_colormap_curve(QPolygonF &curve, qreal x, int iy) // 0.0<=(x)<=1.0, 0<=(iy)<=255
{
    x = qMax(0.0, qMin(1.0,  x));
    qreal y = qMax(0.0, qMin(1.0,  iy/255.0));
    curve << QPointF(x, y);
}
inline void set_colormap_curve(QPolygonF &curve, qreal x, qreal y) // 0.0<=(x, y)<=1.0
{
    x = qMax(0.0, qMin(1.0,  x));
    y = qMax(0.0, qMin(1.0,  y));
    curve << QPointF(x, y);
}
};

#endif // RENDERER_EMT_H
