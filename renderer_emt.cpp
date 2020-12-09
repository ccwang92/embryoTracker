#include "renderer_emt.h"
#include "emt_glwidget.h"
#include <sstream>
#include <string>
#include <cmath>
#include <QMessageBox>
#include <QObject>
#include <QOpenGLExtraFunctions>
//#include <QOpenGLFunctions>

Renderer_EmT::Renderer_EmT(void *widget) : widget(widget)
{
    qDebug(" Renderer::Renderer >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
    //renderer_gl1
    //this->isTera = false; // added by MK, 2018 May, for arranging segments before entering Rnderer_gla::loopCheck
    //this->isLoadFromFile = false;
    //this->pressedShowSubTree = false;
    this->zThick = 1;
    //this->FragTraceMarkerDetector3Dviewer = false;
    //this->NAeditingMode = false;

    //Renderer_EmT
    //tryTexStream = (supported_PBO()? 1 : 0);
    //tryVolShader = (supported_GLSL()? 1 : 0);


    init_members(); //for renderer, gl1 and gl2
}
void Renderer_EmT::makeCurrent()
{
    if (! widget)  return;

    EmT_GLWidget* w = (EmT_GLWidget*)widget;
//	QGLContext* ctx = (QGLContext*)w->context();
//	if ( ctx && ctx->isValid() )
        w->makeCurrent();
}
void Renderer_EmT::initialize()
{
    qDebug("   Renderer_emt::initialize (%d)");

    try {
        //qDebug("  Renderer_gl1::initialize (%d)", version);
        //if (b_error) return; //080924 try to catch the memory error
        // renderer
        ////////////////////////////////////////////////
        //GLeeInit();

        //if (GLEE_ARB_multisample)
        {
            glEnable(GL_MULTISAMPLE_ARB); // default enabled by setSampleBuffers?
            GLint samples=0;
            glGetIntegerv(GL_SAMPLES_ARB, &samples);
            if (samples==2)  glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST); // 5tap for 2sample, 9tap for 4sample
        }

        // GL_*_SMOOTH affected by GL_MULTISAMPLE.
        // Otherwise combined with glEnable(GL_BLEND) & glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE), that is very limited, by RZC 081001
        // Polygon antialiasing is optimized using a blend function (GL_SRC_ALPHA_SATURATE, GL_ONE) with polygons sorted from nearest to farthest.
        //      Destination alpha bit planes, which must be present for this blend function to operate correctly, store the accumulated coverage.
        //glEnable(GL_LINE_SMOOTH);
        //glEnable(GL_POLYGON_SMOOTH); // internal of quads, no effect?

        glDisable (GL_DITHER);
        glDisable (GL_FOG);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        //glEnable(GL_POLYGON_OFFSET_FILL);
        //glPolygonOffset(0, +1); // deal z-fighting, 081121

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);

        if (rgbaBuf==0)   color_background = color_background2; // only geometric objects, 081023

        loadObj();

        // gl1
        loadVol();
        BoundingBox& sBB =surfBoundingBox;
        BoundingBox& BB  =boundingBox;

        BB.x0 = 0;
        BB.y0 = 0;
        BB.z0 = 0;
        BB.x1 = imageX / sampleScaleX * thicknessX;
        BB.y1 = imageY / sampleScaleY * thicknessY;
        BB.z1 = imageZ / sampleScaleZ * thicknessZ;
        qDebug("	BoundingBox surface (%g %g %g)--(%g %g %g)", sBB.x0,sBB.y0,sBB.z0, sBB.x1,sBB.y1,sBB.z1 );
        qDebug("	BoundingBox default (%g %g %g)--(%g %g %g)", BB.x0,BB.y0,BB.z0, BB.x1,BB.y1,BB.z1 );

        loadShader();

    } CATCH_handler( "Renderer_emt::initialize" );
}
void Renderer_EmT::loadVol()
{
    cleanVol(); // 081006: move to before setting imageX/Y/Z, 090705 move to first line
    cleanTexStreamBuffer(); //091012

    qDebug("  Renderer_gl1::loadVol");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop, 081105

    if (! rgbaBuf || bufSize[3]<1 ) return; // no image data, 081002

    ////////////////////////////////////////////////////////////////
    // set coordinate frame size
    sampleScaleX = sampleScale[0];
    sampleScaleY = sampleScale[1];
    sampleScaleZ = sampleScale[2];
    imageX = bufSize[0];
    imageY = bufSize[1];
    imageZ = bufSize[2];
    imageT = bufSize[4];

    bool ok = 0;
    //if ( !(ok = supported_TexNPT()) )
        tryTexNPT = 0;
        qDebug()<< QString("	ARB_texture_non_power_of_two          %1 supported ").arg(ok?"":"NOT");

    //if ( !(ok = supported_TexCompression()) )
        tryTexCompress = 0;
        qDebug()<< QString("	ARB_texture_compression               %1 supported ").arg(ok?"":"NOT");

    //if ( !(ok = supported_Tex3D()) )
        tryTex3D = 0;
        qDebug()<< QString("	EXT_texture3D (or OpenGL 2.0)         %1 supported ").arg(ok?"":"NOT");

    //if ( !(ok = supported_TexStream()) )
        if (tryTexStream != -1)
            tryTexStream = 0;
        qDebug()<< QString("	texture stream (need PBO and GLSL)    %1 supported ").arg(ok?"":"NOT");

    //ok = supported_GL2();
        qDebug()<< QString("	GLSL (and OpenGL 2.0)                 %1 supported ").arg(ok?"":"NOT");


    if (imageT>1) //090802: TexSubImage conflicts against compressed texture2D, but is good for compressed texture3D
    {
        //tryTex3D = 1; 			qDebug("	Turn on tryTex3D for Time series");
        tryTexCompress = 0;		qDebug("		Turn off tryTexCompress for time series");
        tryTexStream = 0;		qDebug("		Turn off tryTexStream for time series");
    }

//	// comment for easy test on small volume
//	if (IS_FITTED_VOLUME(imageX,imageY,imageZ))
//	{
//		if (tryTexStream==1)
//		{
//			qDebug("	No need texture stream for small volume (fitted in %dx%dx%d)", LIMIT_VOLX,LIMIT_VOLY,LIMIT_VOLZ);
//			tryTexStream = 0;  // no need stream, because volume can be fitted in video memory
//		}
//	}

    ////////////////////////////////////////////////////////////////
    // coordinate system
    //
    //     y_slice[z][x]
    //      |
    //      |
    //      |_______ x_slice[z][y]
    //     /
    //    /z_slice[y][x]
    //
    ////////////////////////////////////////////////////////////////
    QTime qtime;  qtime.start();
    qDebug("   setupStack start --- try %s", try_vol_state());

    fillX = _getTexFillSize(imageX);
    fillY = _getTexFillSize(i#include <QOpenGLExtraFunctions>mageY);
    fillZ = _getTexFillSize(imageZ);
    qDebug("   sampleScale = %gx%gx%g""   sampledImage = %dx%dx%d""   fillTexture = %dx%dx%d",
            sampleScaleX, sampleScaleY, sampleScaleZ,  imageX, imageY, imageZ,  fillX, fillY, fillZ);

    //if (tryTex3D && supported_Tex3D())
    {
            qDebug() << "Renderer_gl1::loadVol() - creating 3D texture ID\n";
        glGenTextures(1, &tex3D);		//qDebug("	tex3D = %u", tex3D);
    }
//    if (!tex3D || tryTexStream !=0) //stream = -1/1/2
//    {
//        //tryTex3D = 0; //091015: no need, because tex3D & tex_stream_buffer is not related now.

//            qDebug() << "Renderer_gl1::loadVol() - creating data structures for managing 2D texture slice set\n";

//        Ztex_list = new GLuint[imageZ+1]; //+1 for pbo tex
//        Ytex_list = new GLuint[imageY+1];
//        Xtex_list = new GLuint[imageX+1];
//        memset(Ztex_list, 0, sizeof(GLuint)*(imageZ+1));
//        memset(Ytex_list, 0, sizeof(GLuint)*(imageY+1));
//        memset(Xtex_list, 0, sizeof(GLuint)*(imageX+1));
//        glGenTextures(imageZ+1, Ztex_list);
//        glGenTextures(imageY+1, Ytex_list);
//        glGenTextures(imageX+1, Xtex_list);

//        CHECK_GLErrorString_throw(); // can throw const char* exception, RZC 080925

//        int X = _getBufFillSize(imageX);
//        int Y = _getBufFillSize(imageY);
//        int Z = _getBufFillSize(imageZ);
//        Zslice_data = new RGBA8 [Y * X];//[Z][y][x] //base order
//        Yslice_data = new RGBA8 [Z * X];//[Y][z][x]
//        Xslice_data = new RGBA8 [Z * Y];//[X][z][y]
//        memset(Zslice_data, 0, sizeof(RGBA8)* (Y * X));
//        memset(Yslice_data, 0, sizeof(RGBA8)* (Z * X));
//        memset(Xslice_data, 0, sizeof(RGBA8)* (Z * Y));

//        // optimized copy slice data in setupStackTexture, by RZC 2008-10-04
//    }

    qDebug("   setupStack: id & buffer ....................... cost time = %g sec", qtime.elapsed()*0.001);


    ///////////////////////////////////////
//    if (texture_format==-1)
//    {
//        texture_format = GL_RGBA;
//        //Call TexImage with a generic compressed internal format. The texture image will be compressed by the GL, if possible.
//        //Call CompressedTexImage to Load pre-compressed image.
//        //S3TC: DXT1(1bit alpha), DXT3(sharp alpha), DXT5(smooth alpha)
//        //glHint(GL_TEXTURE_COMPRESSION_HINT_ARB, GL_NICEST); // seems no use, choice DXT3, but DXT5 is better, 081020
//        if (tryTexCompress && GLEE_ARB_texture_compression)
//            texture_format = GL_COMPRESSED_RGBA_ARB;
//        if (texture_format==GL_COMPRESSED_RGBA_ARB && GLEE_EXT_texture_compression_s3tc)
//            texture_format = GL_COMPRESSED_RGBA_S3TC_DXT5_EXT;
//    }
    if (image_format==-1)
    {
        image_format = GL_RGBA;
    }
    if (image_type==-1)
    {
        image_type = GL_UNSIGNED_BYTE;
    }

    subloadTex(volTimePoint, true);   // first loading
    ///////////////////////////////////////

    //update control panel (no need)
    //EmT_GLWidget* w = (EmT_GLWidget*)widget;
    //if (w)  w->updateControl();
}
void Renderer_EmT::subloadTex(V3DLONG timepoint, bool bfirst)
{
    if (texture_format==-1)  return; // not done by loadVol
    if (! rgbaBuf || bufSize[3]<1 ) return; // no image data, 081002

    QTime qtime;  qtime.start();
    {
        timepoint = CLAMP(0, imageT-1, timepoint);
        rgbaBuf = total_rgbaBuf + timepoint*(imageZ*imageY*imageX);

          qDebug() << "Calling setupStackTexture() from Renderer_gl1::subloadTex()";
        //if (tryTexStream<=0) 			// 091014: mix down-sampled & streamed method
        //    setupStackTexture(bfirst);  // use a temporary buffer, so first

        if (tryTexStream >0 && bfirst)
        {
            setupTexStreamBuffer();
        }

    }
    qDebug("   subloadTex [%d]: %s ...... cost time = %g sec", timepoint, try_vol_state(), qtime.elapsed()*0.001);
    qDebug("	  tex_stream_buffer = %s", tex_stream_buffer?"true":"false");

}
int  Renderer_EmT::_getBufFillSize(int w)
{
    return power_of_two_ceil(w); //always use power_of_two
}
int  Renderer_EmT::_getTexFillSize(int w)
{
    if (tryTexNPT)  return w;
    else            return power_of_two_ceil(w);
}
void Renderer_EmT::cleanVol()
{
    qDebug("   Renderer_gl1::cleanVol");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop, 081105

    if (Zslice_data) delete[] Zslice_data;	Zslice_data = 0;
    if (Yslice_data) delete[] Yslice_data;	Yslice_data = 0;
    if (Xslice_data) delete[] Xslice_data;	Xslice_data = 0;
    if (Fslice_data) delete[] Fslice_data;	Fslice_data = 0;

    // solved by makeCurrent: explicit delete of texture-ids causes other view's texture-ids may be deleted. 081021
    if (tex3D) {
        glDeleteTextures(1, &tex3D);
        tex3D = 0;
    }
    if (Ztex_list) { //+1 for pbo tex
        glDeleteTextures(imageZ+1, Ztex_list);
        delete[] Ztex_list;	Ztex_list = 0;
    }
    if (Ytex_list) {
        glDeleteTextures(imageY+1, Ytex_list);
        delete[] Ytex_list;	Ytex_list = 0;
    }
    if (Xtex_list) {
        glDeleteTextures(imageX+1, Xtex_list);
        delete[] Xtex_list;	Xtex_list = 0;
    }
    if (texFslice) {
        glDeleteTextures(1, &texFslice);
        texFslice = 0;
    }

    texture_format = image_format = image_type = -1;
    imageX = imageY = imageZ = imageT = 0;
    safeX = safeY = safeZ = 0;
    fillX = fillY = fillZ = fillF = 0;
    realX = realY = realZ = realF = 0;
}
void Renderer_EmT::cleanData()
{
    qDebug("   Renderer_gl1::cleanData");
    for (int i=0; i<5; i++)
    {
        sampleScale[i]=1; bufSize[i]=0;
    }

    DELETE_AND_ZERO(total_rgbaBuf);
    rgbaBuf = 0;
    DELETE_AND_ZERO(rgbaBuf_Yzx);
    DELETE_AND_ZERO(rgbaBuf_Xzy);
}

void Renderer_EmT::setupData(void* idep)
{
    cleanData(); //090705
    qDebug("  Renderer_gl1::setupData");

    //PROGRESS_DIALOG("", widget);

    this->_idep = idep;

     //b_limitedsize = (tryTexStream==0); //091022, 100720: down-sampling only if tryTexStream==0

     // creating data for 3dviewer when needed
     bool bLocal = false;

     //if (b_limitedsize)
     //{
     //    qDebug("	Down-sampling to 512x512x256 ");
     //}

     try
     {
         Image4DSimple* image4d = ((DataImporter *)_idep)->image4d;
         if (image4d && image4d->getCDim()>0)
         {
             //bLocal = ((iDrawExternalParameter*)_idep)->b_local;
             //bLimited = ((iDrawExternalParameter*)_idep)->b_use_512x512x256; //091015: no need this, because always can use stream texture

             data_unitbytes = image4d->getUnitBytes();
             data4dp = image4d->getRawData();
             //data4d_uint8 = image4d->data4d_uint8;

             size1=dim1 = image4d->getXDim();
             size2=dim2 = image4d->getYDim();
             size3=dim3 = image4d->getZDim();
             size4=dim4 = image4d->getCDim();
             size5=dim5 = 1;
             if (image4d->getTDim()>1 && image4d->getTimePackType()==TIME_PACK_C)
             {
                 //MESSAGE_ASSERT(image4d->getCDim() >= image4d->getTDim());

                 size4=dim4 = image4d->getCDim()/image4d->getTDim();
                 size5=dim5 = image4d->getTDim();
             }
             start1 = 0;
             start2 = 0;
             start3 = 0;
             start4 = 0;
             start5 = 0;

//             if (bLocal)
//             {
//                 size1 = ((iDrawExternalParameter*)_idep)->local_size.x;
//                 size2 = ((iDrawExternalParameter*)_idep)->local_size.y;
//                 size3 = ((iDrawExternalParameter*)_idep)->local_size.z;
//                 sampleScale[0] = sampleScale[1] =sampleScale[2] = 1;

//                 start1 = ((iDrawExternalParameter*)_idep)->local_start.x;
//                 start2 = ((iDrawExternalParameter*)_idep)->local_start.y;
//                 start3 = ((iDrawExternalParameter*)_idep)->local_start.z;
//                 //data4dp += start3*(dim2*dim1) + start2*(dim1) + start1;
//             }
         }
         else // image4d==0  default coordinate frame for surface
         {
             size1=dim1 = 0; //DEFAULT_DIM1;
             size2=dim2 = 0; //DEFAULT_DIM2;
             size3=dim3 = 0; //DEFAULT_DIM3;
             size4=dim4 = 0; // this make no rgbaBuf allocated
             size5=dim5 = 0; // packed time
             start1 = 0;
             start2 = 0;
             start3 = 0;
             start4 = 0;
             start5 = 0;
             data4dp = 0; // this make no rgbaBuf allocated
         }


         bufSize[0] = size1;
         bufSize[1] = size2;
         bufSize[2] = size3;
         bufSize[3] = size4;
         bufSize[4] = size5;
         sampleScale[0]=sampleScale[1]=sampleScale[2]=sampleScale[3]=sampleScale[4] = 1;


         total_rgbaBuf = rgbaBuf = 0; //(RGBA*)-1; //test whether the new sets pointer to 0 when failed
         if (data4dp && size4>0)
         {
             // only RGB, first 3 channels of original image
             total_rgbaBuf = rgbaBuf = new RGBA8[ bufSize[0] * bufSize[1] * bufSize[2] * 1 * bufSize[4] ];
         }

         qDebug("   data4dp = 0x%p \t(start %ldx%ldx%ld_%ld_%ld, size %ldx%ldx%ld_%ld_%ld)", data4dp,
                start1,start2,start3,start4,start5,  size1,size2,size3,size4, size5);

         qDebug("   rgbaBuf = 0x%p \t(%ldx%ldx%ld_%ld_%ld)", rgbaBuf, bufSize[0],bufSize[1],bufSize[2],bufSize[3],bufSize[4]);


         dataViewProcBox = dataBox = BoundingBox(start1, start2, start3, start1+(size1-1), start2+(size2-1), start3+(size3-1));

         qDebug("   data box in original image space @\t(%g %g %g)--(%g %g %g)", dataBox.x0,dataBox.y0,dataBox.z0, dataBox.x1,dataBox.y1,dataBox.z1);

     } CATCH_handler( "Renderer_gl1::setupData" );


     QTime qtime;  qtime.start();
     {
         Image4DSimple* image4d = ((DataImporter *)_idep)->image4d;
         if (image4d)
         {
             Image4DProxy<Image4DSimple> img4dp( image4d );
             img4dp.set_minmax(((DataImporter *)_idep)->p_vmin, ((DataImporter *)_idep)->p_vmax);

             data4dp_to_rgba3d(img4dp,  dim5,
                               start1, start2, start3, start4,
                               size1, size2, size3, size4,
                               total_rgbaBuf, bufSize);
         }

         if (dim4==1)   rgba3d_r2gray(total_rgbaBuf, bufSize); //081103
     }
     qDebug("   data4dp_to_rgba3d ............................................... cost time = %g sec", qtime.elapsed()*0.001);
     // end creating data for 3dviewer
}
QString resourceTextFile(QString filename)
{
    //QFile inputFile(":/subdir/input.txt");
    qDebug() << "Load shader: " << filename;

    QFile inputFile(filename);
    if (inputFile.open(QIODevice::ReadOnly)==false)
        qDebug() << "   *** ERROR in Load shader: " << filename;

    QTextStream in(&inputFile);
    QString line = in.readAll();
    inputFile.close();
    return line;
}
static void linkGLShader(cwc::glShaderManager& SMgr,
        cwc::glShader*& shader, //output
        const char* vertex, const char* fragment)
{
    glGetError(); // clear error
    shader = SMgr.loadfromMemory(vertex, fragment);
    if (shader==0)
    {
       qDebug() << "Renderer_EmT::init:  Error Loading, compiling or linking shader\n";
       throw SMgr.getInfoLog();
    }
}
void Renderer_EmT::loadShader()
{
    cleanShader(); //090705
    qDebug("   Renderer_EmT::loadShader");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop

    try {

        qDebug("+++++++++ shader for Surface Object");
        linkGLShader(SMgr, shaderObj, //0,0);
                //0,
                Q_CSTR(resourceTextFile(":/shader/lighting.txt") + resourceTextFile(":/shader/color_vertex.txt")),
                //Q_CSTR(resourceTextFile(":/shader/vertex_normal.txt")),
                //0);
                Q_CSTR(resourceTextFile(":/shader/lighting.txt") + resourceTextFile(":/shader/obj_fragment.txt")));

        #ifdef TEX_LOD
            QString deftexlod = "#define TEX_LOD \n";
        #else
            QString deftexlod = "#undef TEX_LOD \n";
        #endif

        qDebug("+++++++++ shader for Volume texture2D");
        linkGLShader(SMgr, shaderTex2D,
                0, //Q_CSTR(resourceTextFile(":/shader/color_vertex.txt")),
                Q_CSTR(QString("#undef TEX3D \n") + deftexlod + resourceTextFile(":/shader/tex_fragment.txt")));

        qDebug("+++++++++ shader for Volume texture3D");
        linkGLShader(SMgr, shaderTex3D,
                0, //Q_CSTR(resourceTextFile(":/shader/color_vertex.txt")),
                Q_CSTR(QString("#define TEX3D \n") + deftexlod + resourceTextFile(":/shader/tex_fragment.txt")));

    } CATCH_handler("Renderer_EmT::initialze");
    qDebug("+++++++++ GLSL shader setup finished.");


    glGenTextures(1, &texColormap);
    initColormap();
}

void Renderer_EmT::cleanShader()
{
    qDebug("    Renderer_EmT::cleanShader");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop

    DELETE_AND_ZERO(shaderTex2D);
    DELETE_AND_ZERO(shaderTex3D);
    DELETE_AND_ZERO(shaderObj);

    if (texColormap) {
        glDeleteTextures(1, &texColormap);
        texColormap = 0;
    }
}

void Renderer_EmT::toggleShader()
{
    bool b = (tryVolShader || tryObjShader);
    tryVolShader = !b;
    tryObjShader = !b;
    if (! supported_GLSL())
    {
        qDebug( "	No GL shading language support");
        tryVolShader = 0;
        tryObjShader = 0;
    }
    //qDebug( "	tryShader(vol obj) = %d %d", tryVolShader, tryObjShader);
}
void Renderer_EmT::toggleObjShader()
{
    tryObjShader = !tryObjShader;
    if (! supported_GLSL())
    {
        qDebug( "	No GL shading language support");
        tryObjShader = 0;
    }
    //qDebug( "	tryShader(vol obj) = %d %d", tryVolShader, tryObjShader);
}

void Renderer_EmT::shaderTexBegin(bool stream)
{
    shader = (texture_unit0_3D && !stream)? shaderTex3D : shaderTex2D;

    int format_bgra = (stream && pbo_image_format==GL_BGRA)? 1:0;

    if (IF_VOL_SHADER)
    {
        shader->begin(); //must before setUniform
        shader->setUniform1i("volume",   0); //GL_TEXTURE0
        shader->setUniform1i("colormap", 1); //GL_TEXTURE1

        float n = FILL_CHANNEL-1; // 0-based
        shader->setUniform3f("channel", 0/n, 1/n, 2/n);
        shader->setUniform1i("blend_mode", renderMode);
        shader->setUniform1i("format_bgra", format_bgra);

        // switch to colormap texture
        glActiveTextureARB(GL_TEXTURE1_ARB);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, texColormap);
        glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); // GLSL will replace TexEnv
        CHECK_GLError_print();

//		glTexImage2D(GL_TEXTURE_2D, // target
//				0, // level
//				GL_RGBA, // texture format
//				256, // width
//				FILL_CHANNEL,   // height
//				0, // border
//				GL_RGBA, // image format
//				GL_UNSIGNED_BYTE, // image type
//				&colormap[0][0]);
        glTexSubImage2D(GL_TEXTURE_2D, // target
                0, // level
                0,0, // offset
                256, // width
                FILL_CHANNEL,   // height
                GL_RGBA, // image format
                GL_UNSIGNED_BYTE, // image type
                &colormap[0][0]);
        CHECK_GLError_print();

        // switch back to volume data texture
        glActiveTextureARB(GL_TEXTURE0_ARB);
    }
}

void Renderer_EmT::shaderTexEnd()
{
    if (IF_VOL_SHADER)
    {
        // off colormap
        glActiveTextureARB(GL_TEXTURE1_ARB);
        glDisable(GL_TEXTURE_2D);
        glActiveTextureARB(GL_TEXTURE0_ARB);

        shader->end();
    }
    shader = 0;
}

void Renderer_EmT::equAlphaBlendingProjection()
{
    if (! tryVolShader)
    {
        Renderer_gl1::equAlphaBlendingProjection();
        return;
    }
    glBlendEquationEXT(GL_FUNC_ADD_EXT);/////////////////////
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);	// alpha multiplied in shader
}
void Renderer_EmT::equMaxIntensityProjection()
{
    if (! tryVolShader)
    {
        Renderer_gl1::equMaxIntensityProjection();
        return;
    }
    glBlendEquationEXT(GL_MAX_EXT);//////////////////////////
}

void Renderer_EmT::equMinIntensityProjection()
{
    if (! tryVolShader)
    {
        Renderer_gl1::equMinIntensityProjection();
        return;
    }
    glBlendEquationEXT(GL_MIN_EXT);//////////////////////////
}

void Renderer_EmT::equCrossSection()
{
    if (! tryVolShader)
    {
        Renderer_gl1::equCrossSection();
        return;
    }
    glBlendEquationEXT(GL_FUNC_ADD_EXT);/////////////////////
    glBlendColorEXT(1, 1, 1, 1-CSbeta);
    glBlendFunc(GL_CONSTANT_ALPHA, GL_ONE_MINUS_CONSTANT_ALPHA); // constant Alpha
}

void Renderer_EmT::initColormap()
{
    qDebug("   Renderer_EmT::initColormap");

    for (int ch=0; ch<N_CHANNEL; ch++)
    {
        for (int i=0; i<256; i++) // intensity
        {
            for (int j=0; j<3; j++)  colormap[ch][i].c[j] = (ch==j)*255;	colormap[ch][i].a = i;
            //for (int j=0; j<3; j++)  colormap[ch][i].c[j] = (ch!=j)*i;	colormap[ch][i].a = i;
        }

        for (int j=0; j<4; j++) // RGBA
        {
            colormap_curve[ch][j].clear();
            int y;
            y = colormap[ch][0].c[j];	   set_colormap_curve(colormap_curve[ch][j],  0.0,  y);
            y = colormap[ch][255].c[j];    set_colormap_curve(colormap_curve[ch][j],  1.0,  y);

        //	qDebug() << QString("[%1][%2]").arg(ch).arg(j) <<  colormap_curve[ch][j];
        }
    }

    glBindTexture(GL_TEXTURE_2D, texColormap);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); // MUST use nearest filter
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); // MUST use nearest filter

    // load texture
    glTexImage2D(GL_TEXTURE_2D, // target
            0, // level
            GL_RGBA, // texture format
            256, // width
            FILL_CHANNEL,   // height
            0, // border
            GL_RGBA, // image format
            GL_UNSIGNED_BYTE, // image type
            &colormap[0][0]);
            //NULL); ////////////////////////////////////// load ON DEMAND when drawing
    CHECK_GLError_print();

}

void Renderer_EmT::applyColormapToImage() // at most the first 3 channels
{
    qDebug("   Renderer_EmT::applyColormapToImage");

#ifndef test_main_cpp
    My4DImage* curImg = v3dr_getImage4d(_idep);
    if (curImg && data4dp && dim4>0)
    {
        qDebug("	dim[%dx%dx%d_%d]", dim1,dim2,dim3,dim4);
        int ix, iy, iz;
        for (iz = 0; iz < dim3; iz++)
            for (iy = 0; iy < dim2; iy++)
                for (ix = 0; ix < dim1; ix++)
                {
                    RGB8 oldC = getRGB3dUINT8 (data4dp, dim1, dim2, dim3, dim4,  ix, iy, iz);

                    RGB8 newC = lookupColormap(oldC, 1); //OP_ADD

                    setRGB3dUINT8 (data4dp, dim1, dim2, dim3, dim4,  ix, iy, iz,  newC);
                }

        curImg->updateViews();
    }
#endif
}

////////////////////////////////////////////////////
// shader for cross-section mode in tex_fragment.txt
//
//vec3 aC1 = C1.rgb * C1.a;
//vec3 aC2 = C2.rgb * C2.a;
//vec3 aC3 = C3.rgb * C3.a;
//float Amean = (C1.a + C2.a + C3.a)/3.0;
//			// pow((C1.a * C2.a * C3.a), 1.0/3.0);
//float Amax	= max(C1.a, max(C2.a, C3.a));
//
//vec4 oColor;
//if (mode==0) // cross-section
//{
//	float Axsec = Amean;
//	oColor.rgb = (aC1 + aC2 + aC3);
//	oColor.a = Axsec;
//}

#define OP_MAX	0
#define OP_ADD	1

RGB8 Renderer_EmT::lookupColormap(RGB8 inC, int op)
{
    #define R(k,j)	(colormap[k][j].r/255.0)
    #define G(k,j)	(colormap[k][j].g/255.0)
    #define B(k,j)	(colormap[k][j].b/255.0)
    #define A(k,j)	(colormap[k][j].a/255.0)

    #define AR(k,j)	(A(k,j)*R(k,j))
    #define AG(k,j)	(A(k,j)*G(k,j))
    #define AB(k,j)	(A(k,j)*B(k,j))

    int i1 = inC.r;
    int i2 = inC.g;
    int i3 = inC.b;

    float o1,o2,o3; // o1=o2=o3=0;

    if (op==OP_MAX)
    {
        o1 = MAX(AR(1,i1), MAX(AR(2,i2), AR(3,i3)));
        o2 = MAX(AG(1,i1), MAX(AG(2,i2), AG(3,i3)));
        o3 = MAX(AB(1,i1), MAX(AB(2,i2), AB(3,i3)));

    }
    else if (op==OP_ADD)
    {
        o1 = AR(1,i1) + AR(2,i2) + AR(3,i3);
        o2 = AG(1,i1) + AG(2,i2) + AG(3,i3);
        o3 = AB(1,i1) + AB(2,i2) + AB(3,i3);
    }

    RGB8 oC;
    oC.r = o1*255;
    oC.g = o2*255;
    oC.b = o3*255;
    return oC;
}

bool Renderer_EmT::supported_TexStream()
{
    if (imageT>1)
    {
        qDebug( "		Time series is NOT supported by texture stream!");
        tryTexStream = 0;
        return false;
    }
//	if (sizeof(void*)<8)
//	{
//		qDebug( "		32-bit system is NOT supported by texture stream!");
//		tryTexStream = 0;
//		return false;
//	}
    if (!supported_PBO())
    {
        qDebug( "		ARB_pixel_buffer_object 	NOT supported !");
        tryTexStream = 0;
        return false;
    }
    if (!supported_GLSL())
    {
        qDebug( "		ARB_shading_language_100 	NOT supported !");
        tryTexStream = 0;
        return false;
    }
    return true;
}

/////////////////////////////////////////////////////////////////////
//Streaming textures using pixel buffer objects:
//
//    const int texWidth = 256;
//    const int texHeight = 256;
//    const int texsize = texWidth * texHeight * 4;
//    void *pboMemory, *texData;
//
//    // Define texture level zero (without an image); notice the
//    // explicit bind to the zero pixel unpack buffer object so that
//    // pass NULL for the image data leaves the texture image
//    // unspecified.
//    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, texWidth, texHeight, 0,
//                 GL_BGRA, GL_UNSIGNED_BYTE, NULL);
//
//    // Create and bind texture image buffer object
//    glGenBuffers(1, &texBuffer);
//    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, texBuffer);
//
//    // Setup texture environment
//    ...
//
//    texData = getNextImage();
//
//    while (texData) {
//
//        // Reset the contents of the texSize-sized buffer object
//        glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, texSize, NULL,
//                     GL_STREAM_DRAW);
//
//        // Map the texture image buffer (the contents of which
//        // are undefined due to the previous glBufferData)
//        pboMemory = glMapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB,
//                                GL_WRITE_ONLY);
//
//        // Modify (sub-)buffer data
//        memcpy(pboMemory, texData, texsize);
//
//        // Unmap the texture image buffer
//        glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER_ARB);
//
//        // Update (sub-)teximage from texture image buffer
//        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, texWidth, texHeight,
//                        GL_BGRA, GL_UNSIGNED_BYTE, BUFFER_OFFSET(0));
//
//        // Draw textured geometry
//        glBegin(GL_QUADS);
//        ...
//        glEnd();
//
//        texData = getNextImage();
//    }
//
//    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

#define BIND_UNPACK_PBO(pbo)  glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo)


void Renderer_EmT::cleanTexStreamBuffer()
{
    qDebug("    Renderer_EmT::cleanStreamBuffer");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop

    // release PBO object
    BIND_UNPACK_PBO(0);
    if (pboZ) {
        glDeleteBuffersARB(1, &pboZ);
        pboZ = 0;
    }
    if (pboY) {
        glDeleteBuffersARB(1, &pboY);
        pboY = 0;
    }
    if (pboX) {
        glDeleteBuffersARB(1, &pboX);
        pboX = 0;
    }
    tex_stream_buffer = pboZ>0; //#################

    DELETE_AND_ZERO(rgbaBuf_Yzx);
    DELETE_AND_ZERO(rgbaBuf_Xzy);
}

const bool use_pboZ_only = false;

void  Renderer_EmT::setupTexStreamBuffer()
{
    cleanTexStreamBuffer();
    qDebug("   Renderer_EmT::setupStreamBuffer");
    makeCurrent(); //ensure right context when multiple views animation or mouse drop

    //MESSAGE_ASSERT( !tryTexCompress && !tryTex3D && imageT<=1 );
    try
    {
        DELETE_AND_ZERO(rgbaBuf_Yzx);
        DELETE_AND_ZERO(rgbaBuf_Xzy);
    //	rgbaBuf_Yzx = new RGBA8 [imageX*imageY*imageZ];  // this direction copy does not help for speed
    //	_copyYzxFromZyx( rgbaBuf_Yzx, rgbaBuf,  imageX, imageY, imageZ);
        rgbaBuf_Xzy = new RGBA8 [imageX*imageY*imageZ];
        _copyXzyFromZyx( rgbaBuf_Xzy, rgbaBuf,  imageX, imageY, imageZ);
    }
    catch(...)
    {
        qDebug("	setupTexStreamBuffer: Out of memory. It may cause rendering speed down!");
    }
    fillX = _getTexFillSize(imageX);
    fillY = _getTexFillSize(imageY);
    fillZ = _getTexFillSize(imageZ);


    //091012: 1 common PBO is not faster than switch 3 PBOs
    BIND_UNPACK_PBO(0);
    glGenBuffers(1, &pboZ);
    glGenBuffers(1, &pboY);
    glGenBuffers(1, &pboX);
    tex_stream_buffer = pboZ>0; //##################

    glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      //// 4-byte pixel alignment image for good speed
    pbo_texture_format = GL_RGBA8;
    //pbo_texture_format = GL_RGB8;				//// no help for speed
    //pbo_image_format = GL_RGBA;
    pbo_image_format = GL_BGRA;  					//// only BGRA format is speeded up by PBO
    pbo_image_type = GL_UNSIGNED_BYTE;				//// no care about big-endian or little-endian
    //pbo_image_type = GL_UNSIGNED_INT_8_8_8_8_REV;
    // When using GL_UNSIGNED_INT_8_8_8_8_REV, the OpenGL implementation expects to find data in byte order ARGB on big-endian systems, but BGRA on little-endian systems.
    // for GeForce 8800 GT 512M PCIE x16, 1GB RGB8 image, speed is 0.3hz(no compression, no pbo), 0.6hz(pbo RGBA), 1.1hz(pbo BGRA), 2.2hz(compression, no pbo)

    int max_w = MAX(fillX, MAX(fillY,fillZ));
    for (int stack_i=1; stack_i<=3; stack_i++)
    {
        GLuint pbo = 0, tex = 0;
        int w = 0, h =0;
        int size = 0;

        switch (stack_i)
        {
        case 1: //Z[y][x]
            pbo = pboZ;
            tex = Ztex_list[0];
            w = fillX, h = fillY;
            break;
        case 2: //Y[z][x]
            pbo = pboY;
            tex = Ytex_list[0];
            w = fillX, h = fillZ;
            break;
        case 3: //X[z][y]
            pbo = pboX;
            tex = Xtex_list[0];
            w = fillY, h = fillZ;
            break;
        }

        //MUST pre-allocated PBO memory size not less than used, otherwise report invalid operation
        size = w*h*4;
        if (use_pboZ_only)
        {
            size = max_w*max_w*4;
            pbo = pboZ; // use 1 common PBO
        }
        ///////////////////////////////////////////
        BIND_UNPACK_PBO(pbo);
        glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, size, NULL, GL_STREAM_DRAW);

        glBindTexture(GL_TEXTURE_2D, tex);
        setTexParam2D();
        glTexImage2D(GL_TEXTURE_2D, // target
            0, // level
            pbo_texture_format, // texture format
            w, // width
            h, // height
            0, // border
            pbo_image_format, // image format
            pbo_image_type, // image type
            NULL);
        CHECK_GLErrorString_throw(); // can throw const char* exception, RZC 080925
    }
    BIND_UNPACK_PBO(0);
}

void Renderer_EmT::setupStackTexture(bool bfirst) {
    qDebug() << "Renderer_EmT::setupStackTexture() start";
    Renderer_EmT::setupStackTexture(bfirst);
}


void Renderer_EmT::_streamTex(int stack_i, int slice_i, int step, int slice0, int slice1)
{
    GLuint pbo = 0, tex = 0;
    RGBA8 *pbo_mem = 0, *p_slice = 0;
    int w = 0, h =0;
    int sw = 0, sh = 0;
    int size = 0;

    switch (stack_i)
    {
    case 1: //Z[y][x]
        pbo = pboZ;
        tex = Ztex_list[0];
        p_slice = Zslice_data;
        w = fillX, h = fillY;
        sw = COPY_X, sh = COPY_Y;
        break;
    case 2: //Y[z][x]
        pbo = pboY;
        tex = Ytex_list[0];
        p_slice = Yslice_data;
        w = fillX, h = fillZ;
        sw = COPY_X, sh = COPY_Z;
        break;
    case 3: //X[z][y]
        pbo = pboX;
        tex = Xtex_list[0];
        p_slice = Xslice_data;
        w = fillY, h = fillZ;
        sw = COPY_Y, sh = COPY_Z;
        break;
    }

    //sw = w; sh =h;
    size = sw*sh*4;
    if (use_pboZ_only)
    {
        pbo = pboZ; // use 1 common PBO
    }
    if (!pbo || !tex) return;

    MESSAGE_ASSERT(imageX>=realX && imageY>=realY && imageZ>=realZ);
    MESSAGE_ASSERT(COPY_X>=realX && COPY_Y>=realY && COPY_Z>=realZ);
//	sw+=1, sh+=1;  // to get rid of artifacts
//	if (sw>w) sw = w;
//	if (sh>h) sh = h;

    glBindTexture(GL_TEXTURE_2D, tex);
    BIND_UNPACK_PBO(pbo);
    glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, size, NULL, GL_STREAM_DRAW);
    CHECK_GLError_print();

    pbo_mem = (RGBA8*)glMapBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, GL_WRITE_ONLY);
    CHECK_GLError_print();
    if (pbo_mem)
    {
        //memset(pbo_mem, slice_i, size);
        //memcpy(pbo_mem, rgbaBuf + slice_i, size);
        _copySliceFromStack(rgbaBuf, imageX,imageY,imageZ,  pbo_mem, sw,  stack_i, slice_i,  rgbaBuf_Yzx, rgbaBuf_Xzy);

        // Unmap the texture image buffer & Start DMA transfer
        glUnmapBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB);
        CHECK_GLError_print();
    }
    else
    {
        //memset(p_slice, 0, size);
        _copySliceFromStack(rgbaBuf, imageX,imageY,imageZ,  p_slice, sw,  stack_i, slice_i,  rgbaBuf_Yzx, rgbaBuf_Xzy);
    }

    setTexParam2D(); //100809
    glTexSubImage2D(GL_TEXTURE_2D, // target
        0, // level
        0,0,  // offset
        sw, // width
        sh, // height
        pbo_image_format, // image format
        pbo_image_type, // image type
        (pbo_mem)? 0    // PBO offset
            : p_slice); // CPU memory
    CHECK_GLError_print();
}

void Renderer_EmT::_streamTex_end()
{
    BIND_UNPACK_PBO(0);
}

bool Renderer_EmT::_streamTex_ready()
{
    // tex_stream_buffer is set only when createStreamBuffer is successful.
    // tryTexStream: -1--resident, 0--off, 1--mixed, 2--force
    // beStill: means no user input a while
    bool need_stream = (tryTexStream==2)
                    || (tryTexStream==1 && beStill())
                    || (renderMode==rmCrossSection); //20110709
    return  (tex_stream_buffer && need_stream);
}

void Renderer_EmT::toggleTexStream()
{
    tryTexStream = !(tryTexStream >0); // -1--resident, 0--off, 1--mixed, 2--force
    //qDebug( "	tryTexStream = %d", tryTexStream);
    try	{
        //PROGRESS_DIALOG( ((tryTexStream >0)? "Try Texture Stream": "No Texture Stream"), widget);
        //PROGRESS_PERCENT(30);

        if (tryTexStream<=0)
        {
            cleanTexStreamBuffer();
            qDebug("	toggleTexStream: %s", try_vol_state());
        }
        else loadVol();

        //PROGRESS_PERCENT(100);
    } CATCH_handler( "Renderer_EmT::toggleTexStream" );
}



void Renderer_EmT::setupView(int width, int height)
{
    //qDebug(" Renderer::setupView");
    makeCurrent(); // Qt seems not makeCurrent in resizeGL, 081029 by RZC

//	if (screenW != width || screenH != height)
//	{
//		if (depth_buffer)	delete[] depth_buffer;	depth_buffer = 0;
//		depth_buffer = new GLfloat[width*height];
//	}

    screenW = width;
    screenH = height;
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    setProjection();

    glMatrixMode(GL_MODELVIEW);
    // here please make object space fit in [-1, +1]^3
}
void Renderer_EmT::setProjection()
{
    double aspect = double(screenW)/MAX(screenH,1);

    if (bOrthoView)
    {
        double halfw = 1.3*aspect *zoomRatio;
        double halfh = 1.3        *zoomRatio;
        glOrtho(-halfw, halfw, -halfh, halfh, viewNear, viewFar);
    }
    else
    {
        gluPerspective(viewAngle*zoomRatio, aspect, viewNear, viewFar);
    }

    glTranslated(0, 0, -viewDistance);
}
void Renderer_EmT::paint()
{
    //qDebug(" Renderer_gl1::paint(renderMode=%i)", renderMode);

    if (b_error) return; //080924 try to catch the memory error

    glClearColor(color_background.r, color_background.g, color_background.b, 0);
    glDepthRange(0, 1);

    // CQB 2015/12/16: performance optimization: at high resolutions, drawing in the track lags
    // significantly because this paint routine spends most of its time redrawing the volume.
    // since the volume can't change while the user is drawing a curve, we skip all of the rendering
    // steps except drawing the track while the track is being displayed.
//    if (!sShowTrack || highlightedEndNodeChanged)
//    {
//        glClearStencil(0);
//        glClearDepth(1);
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

//        // clearing framebuffer; reset the drawn flag on all the markers
//        for (int marker = 0; marker < listMarkerPos.size(); marker++)
//        {
//            listMarkerPos[marker].drawn = false;
//        }
//    }

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    //GL_LEQUAL);

    glMatrixMode(GL_MODELVIEW);
    // here, MUST BE normalized space of [-1,+1]^3;

    glGetDoublev(GL_MODELVIEW_MATRIX, volumeViewMatrix); //no scale here, used for drawUnitVolume()

    glPushMatrix();
    //setMarkerSpace(); // space to define marker & curve
    //glGetIntegerv(GL_VIEWPORT,         viewport);            // used for selectObj(smMarkerCreate)
    //glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);    // used for selectObj(smMarkerCreate)
    //glGetDoublev(GL_MODELVIEW_MATRIX,  markerViewMatrix);    // used for selectObj(smMarkerCreate)
    glPopMatrix();

    bShowCSline = bShowAxes;
    bShowFSline = bShowBoundingBox;

    prepareVol();

    if (!sShowTrack)
    {
        if (!b_renderTextureLast) {
            renderVol();
        }

//        if (sShowMarkers>0 || sShowSurfObjects>0)
//        {
//            if (polygonMode==1)	      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//            else if (polygonMode==2)  glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
//            else                      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//            setObjLighting();

//            if (sShowSurfObjects>0)
//            {
//                glPushMatrix(); //================================================= SurfObject {

//                // original surface object space ==>fit in [-1,+1]^3
//                setSurfaceStretchSpace();
//                glPushName(dcSurface);
//                drawObj();  // neuron-swc, cell-apo, label-surf, etc
//                glPopName();

//                glPopMatrix(); //============================================================= }
//            }

//            if (sShowMarkers>0)
//            {
//                glPushMatrix(); //===================================================== Marker {

//                // marker defined in original image space ==>fit in [-1,+1]^3
//                setMarkerSpace();
//                glPushName(dcSurface);
//                drawMarker();  // just markers
//                glPopName();

//                glPopMatrix(); //============================================================= }
//            }

//            disObjLighting();
//        }

//        if (! b_selecting)
//        {
//            if (bOrthoView)
//            {
//                glPushMatrix(); //============================================== scale bar {
//                //drawScaleBar();
//                glPopMatrix(); //========================================================= }
//            }
//        }

        if (b_renderTextureLast) {
            renderVol();
        }

        // must be at last
        // show rubber band track for dragging neuron
//        if (! b_selecting && sShowRubberBand)
//        {
//            if (polygonMode==1)	      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//            else if (polygonMode==2)  glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
//            else                      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//            setObjLighting();
//            beginHighlight();
//            glPushMatrix();
//            //setMarkerSpace();
//            //blendRubberNeuron();

////#ifndef test_main_cpp //140211
////            blendDraggedNeuron();
////#endif

//            glPopMatrix();
//            endHighlight();
//            disObjLighting();
//        }
    } // !sShowTrack

//    if (! b_selecting && sShowTrack)
//    {
//        blendTrack();
//    }

    //always draw some background text by PHC 20151117

    if (1)
    {
        //draw in the center of the display, - disabled for now as it might be disturbing for viewing. PHC 20151117
//        float D = (BB.Dmax());
//        XYZ A0 = BB.Vabsmin();
//        glColor3f(0.7, 0.7, 0.7);
        //drawString(A0.x+BB.Dx()/2, A0.y+BB.Dy()/2, A0.z+BB.Dz()/2, "vaa3d.org", 0, 50);

        //draw at the corner
        glPushMatrix(); //============================================== {

        //drawVaa3DInfo(16);
        //drawEditInfo();


        //drawSegInfo();

        glPopMatrix(); //========================================================= }
    }

    return;
}
void Renderer_EmT::renderVol()
{
    if (has_image())
    {
            glPushMatrix(); //===================================================== Volume {

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

            //if (renderMode==rmAlphaBlendingProjection || renderMode==rmMaxIntensityProjection || renderMode==rmMinIntensityProjection) // not for rmCrossSection
                    enableViewClipPlane(); //front-cut-plane

            // unit image space ==>fit in [-1,+1]^3
            setUnitVolumeSpace();
            glPushName(dcVolume);

            drawVol();
            glPopName();

//            if (! b_selecting) if (bShowCSline && renderMode==rmCrossSection)
//            {
//                    drawCrossLine(2);
//            }
//            if (! b_selecting) if (bShowFSline && bFSlice) //renderMode==rmCrossSection)
//            {
//                    drawUnitFrontSlice(1); // just draw bound Line of F-Slice
//            }

            //if (renderMode==rmAlphaBlendingProjection || renderMode==rmMaxIntensityProjection || renderMode==rmMinIntensityProjection)
                    disableViewClipPlane();

            glPopMatrix(); //============================================================== }
        }

}
void Renderer_EmT::drawVol()
{
    if (volTimeOffset)
    {
        qDebug(" volTimeOffset = %d, volTimePoint = %d", volTimeOffset, volTimePoint);

        volTimePoint = CLAMP(0, imageT-1, volTimePoint);
        volTimeOffset = 0;
        subloadTex(volTimePoint, false);
    }

    //---------------------------------------
    if (CSbeta >=1) return; // 081202, avoid be hit when invisible

    float af = 3.0f;
                //1.0f;
    alpha_threshold = pow(double(CSbeta), double(af));
    color_proxy.a = pow(1.0-double(CSbeta), double(1.0/af)); //change to double() by PHC, 2010-05-20
    SLICE_COLOR = color_proxy;

    if (has_image() && !b_renderTextureLast) // if rendering texture first, we can clear - otherwise this is done in prepareVol()
    {
        glColor3f(0, 0, 0);
        drawBackFillVolCube(); // clear the project region to zero for MIP
    }
    glEnable(GL_BLEND);      equMaxIntensityProjection();
    glEnable(GL_ALPHA_TEST); glAlphaFunc(GL_GEQUAL, alpha_threshold); // >= threshold Alpha, 080930

//	switch (renderMode)
//	{
//	case rmAlphaBlendingProjection:
//		glEnable(GL_BLEND);      equAlphaBlendingProjection();
//		glEnable(GL_ALPHA_TEST); glAlphaFunc(GL_GREATER, alpha_threshold); // > threshold Alpha
//		break;

//	case rmMaxIntensityProjection:
//        if (has_image() && !b_renderTextureLast) // if rendering texture first, we can clear - otherwise this is done in prepareVol()
//		{
//			glColor3f(0, 0, 0);
//			drawBackFillVolCube(); // clear the project region to zero for MIP
//		}
//		glEnable(GL_BLEND);      equMaxIntensityProjection();
//		glEnable(GL_ALPHA_TEST); glAlphaFunc(GL_GEQUAL, alpha_threshold); // >= threshold Alpha, 080930
//		break;

//	case rmMinIntensityProjection:
//		if (has_image() && !b_renderTextureLast) // if rendering texture first, we can clear - otherwise this is done in prepareVol()
//		{
//			glColor3f(0.8, 0.8, 0.8);
//			drawBackFillVolCube(); // clear the project region to a high gray-level value for MIP (black won't do fine)
//		}
//		glEnable(GL_BLEND);      equMinIntensityProjection();
//		glEnable(GL_ALPHA_TEST); glAlphaFunc(GL_LEQUAL, 1 - alpha_threshold); // >= threshold Alpha, 080930
//		break;

//	case rmCrossSection:
//		if (GLEE_EXT_blend_color)
//		{
//			glEnable(GL_BLEND);  equCrossSection();
//		}
//		else
//		{
//			glDisable(GL_BLEND);
//		}
//		glDisable(GL_ALPHA_TEST);
////		if (CSbeta >=0.5)
////			glDisable(GL_DEPTH_TEST);
////		else
////			glEnable(GL_DEPTH_TEST);
//		break;

//	default: // unknown ?
//		break;
//	}

    // modulate by color_proxy, by RZC 080930
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE); // this is default, but make sure;
    //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glColor4fv(SLICE_COLOR.c);

    glDisable(GL_LIGHTING); //volume lighting is difference from polygon object

    glShadeModel(GL_FLAT); //flat for slice rectangle
    {
        updateVolCutRange();

        glEnable(GL_TEXTURE_3D); //glBindTexture(GL_TEXTURE_3D, 0);
        glEnable(GL_TEXTURE_2D); //glBindTexture(GL_TEXTURE_2D, 0);

        glPushName(renderMode);  //100729 make nameLength>2 to skip the Intel GL ICD bug
        {
            drawUnitVolume(); //100729 select name in drawStackX/Y/Z

            glPushName(vsFslice);  //100729 add select name of vsFslice
                drawUnitFrontSlice(0); // F-Slice
            glPopName();
        }
        glPopName();

        glDisable(GL_TEXTURE_3D);
        glDisable(GL_TEXTURE_2D);

    }
    glShadeModel(GL_SMOOTH);

    glDisable(GL_BLEND); //090429 RZC: no effect to glBlendEquationEXT(GL_MAX_EXT), must set to GL_FUNC_ADD_EXT
    glBlendEquationEXT(GL_FUNC_ADD_EXT); //defined in glew.h
    glDisable(GL_ALPHA_TEST);
}
void Renderer_EmT::prepareVol()
{
    // In the b_renderTextureLast case we need to clear the volume before we draw the markers, not after.
    // Note that in the case where the textures are rendered first, drawVol() will
    // clear the volume if MIP is the mode, so we don't have to do it here.
    if (has_image() && b_renderTextureLast)
    {
        glPushMatrix(); //===================================================== Volume {

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

//        if (renderMode==rmAlphaBlendingProjection || renderMode==rmMaxIntensityProjection || renderMode==rmMinIntensityProjection) // not for rmCrossSection
//            enableViewClipPlane(); //front-cut-plane

//        // unit image space ==>fit in [-1,+1]^3
//        setUnitVolumeSpace();

//        if (renderMode==rmMaxIntensityProjection) {
//            glColor3f(0, 0, 0);
//            drawBackFillVolCube(); // clear the project region to zero for MIP
//        }

//        if (renderMode==rmMinIntensityProjection) {
//            glColor3f(0.8, 0.8, 0.8);
//            drawBackFillVolCube(); // clear the project region to near-white for mIP
//        }

        enableViewClipPlane(); //front-cut-plane
        setUnitVolumeSpace();
        glColor3f(0, 0, 0);
        drawBackFillVolCube(); // clear the project region to zero for MIP
        glPopMatrix(); //============================================================== }
    }
}
void Renderer_EmT::enableViewClipPlane()
{
    //qDebug("	 enableViewClipPlane  %g", viewClip);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // in View space
    // keep for dot(clip, pos)>=0
    GLdouble clip0[] = { 0.0,  0.0, -1.0,  0 };
    // [0, 1] ==> [+1, -1]*(s)
    clip0[3] = viewClip;

    //configure the clip planes
    glClipPlane(GL_CLIP_PLANE0, clip0);
    glEnable(GL_CLIP_PLANE0);

    glPopMatrix();
}
void Renderer_EmT::setUnitVolumeSpace()
{
    BoundingBox BB = boundingBox;  // a copy
    float DX = BB.Dx();
    float DY = BB.Dy();
    float DZ = BB.Dz();
    float maxD = BB.Dmax();

    float s[3];
    s[0] = DX /maxD *2;
    s[1] = DY /maxD *2;
    s[2] = DZ /maxD *2;

    // form unit volume space ==> fit in [-1, +1]^3
    glScaled(s[0], s[1], s[2]);
    //qDebug("Scale from [0,1]: x=%f, y=%f, z=%f", s[0],s[1],s[2]);
    glTranslated(-.5, -.5, -.5);
}
void Renderer_EmT::drawBackFillVolCube()
{
    if ((VOL_X1-VOL_X0<0)||(VOL_Y1-VOL_Y0<0)||(VOL_Z1-VOL_Z0<0)) return;

    glPushAttrib(GL_DEPTH_BUFFER_BIT);
    //glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glBegin(GL_QUADS);
    {
        //yx0
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z0);
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z0);
        //x0z
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z1);
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z1);
        //0zy
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z0);
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z1);
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z1);
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z0);
        //xy1
        glVertex3f(VOL_X0, VOL_Y0, VOL_Z1);
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z1);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z1);
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z1);
        //z1x
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z0);
        glVertex3f(VOL_X0, VOL_Y1, VOL_Z1);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z1);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z0);
        //1yz
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z0);
        glVertex3f(VOL_X1, VOL_Y1, VOL_Z1);
        glVertex3f(VOL_X1, VOL_Y0, VOL_Z1);
    }
    glEnd();
    glPopAttrib();
}

void Renderer_EmT::blendBrighten(float fbright, float fcontrast) // fast, 8-bit precision
{
    //fcontrast = 1.5;
    if (fbright==0 && fcontrast==1) return;

    if (fbright>=-1 && fbright<=1)
    {
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_DEPTH_TEST);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        glEnable(GL_BLEND);
#define RECT_BLEND	glRecti(-100,-100, 100,100) // assume [-1,+1]^3 view

        if (fcontrast <=1)
        {
            if (fbright >=0)
            {
                //glBlendEquationEXT(GL_FUNC_ADD_EXT);
                glBlendEquation(GL_FUNC_ADD);
            }
            else //fbright <0
            {
                glBlendEquationEXT(GL_FUNC_REVERSE_SUBTRACT_EXT);
                fbright = -fbright;
            }
            glBlendFunc(GL_ONE, GL_SRC_ALPHA); // new_color = 1*fbright + fcontrast*old_color
            glColor4f(fbright, fbright, fbright, fcontrast);
            RECT_BLEND;
        }

        else //fcontrast >1
        {
            float res = fcontrast -1; // int(fcontrast);
            if (res)
            {
                glBlendEquationEXT(GL_FUNC_ADD_EXT);
                glBlendFunc(GL_DST_COLOR, GL_ONE); // new_color = res*old_color + 1*old_color;
                glColor4f(res, res, res, 1);
                RECT_BLEND;
            }

            if (fbright >=0)
            {
                glBlendEquationEXT(GL_FUNC_ADD_EXT);
            }
            else //fbright <0
            {
                glBlendEquationEXT(GL_FUNC_REVERSE_SUBTRACT_EXT);
                fbright = -fbright;
            }
            glBlendFunc(GL_ONE, GL_ONE); // new_color = 1*fbright + 1*old_color
            glColor4f(fbright, fbright, fbright, 1);
            RECT_BLEND;
        }

        glPopMatrix();
        glPopAttrib();
    }
}


void Renderer_EmT::rgba3d_r2gray(RGBA8* rgbaBuf, V3DLONG bufSize[5])
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

void Renderer_EmT::data4dp_to_rgba3d(Image4DProxy<Image4DSimple>& img4dp, V3DLONG dim5,
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

void Renderer_EmT::data4dp_to_rgba3d(unsigned char* data4dp, V3DLONG dim1, V3DLONG dim2, V3DLONG dim3, V3DLONG dim4, V3DLONG dim5,
        V3DLONG start1, V3DLONG start2, V3DLONG start3, V3DLONG start4,
        V3DLONG size1, V3DLONG size2, V3DLONG size3, V3DLONG size4,
        RGBA8* rgbaBuf, V3DLONG bufSize[5])
{
    if (data4dp==0 || rgbaBuf==0 || bufSize==0)
        return;

    //V3DLONG block_size = dim3*dim2*dim1;
    #define SAMPLE(it, ic, ix,iy,iz, dx,dy,dz) \
                (unsigned char)sampling3dAllTypes( data4dp+ (it*dim4 + ic)*(dim3*dim2*dim1), \
                                                dim1, dim2, dim3, ix, iy, iz, dx, dy, dz )

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
    MESSAGE_ASSERT(dx*dy*dz >=1); //down sampling

    V3DLONG ot;
    V3DLONG ox, oy, oz;
    V3DLONG ix, iy, iz;
    for (ot=0; ot<imageT; ot++)
    for (oz = 0; oz < imageZ; oz++)
    for (oy = 0; oy < imageY; oy++)
    for (ox = 0; ox < imageX; ox++)
        {
            ix = start1+ CLAMP(0,dim1-1, IROUND(ox*sx));
            iy = start2+ CLAMP(0,dim2-1, IROUND(oy*sy));
            iz = start3+ CLAMP(0,dim3-1, IROUND(oz*sz));

            RGBA8 rgba;

            if (imageC >= 1) {
                rgba.r = SAMPLE(ot, 0, ix,iy,iz, dx,dy,dz);
            } else {
                rgba.r = 0;
            }

            if (imageC >= 2) {
                rgba.g = SAMPLE(ot, 1, ix,iy,iz, dx,dy,dz);
            } else {
                rgba.g = 0;
            }

            if (imageC >= 3) {
                rgba.b = SAMPLE(ot, 2, ix,iy,iz, dx,dy,dz);
            } else {
                rgba.b = 0;
            }

            if (imageC >= 4) {
                rgba.a = SAMPLE(ot, 3, ix,iy,iz, dx,dy,dz);
            } else {
                float t = //MAX(rgba.r, MAX(rgba.g, rgba.b));
                            ((0.f + rgba.r + rgba.g + rgba.b) / imageC);
                rgba.a = (unsigned char)t;
                            //(unsigned char)(t*t/255);
                            //(unsigned char)(sqrt(t/255)*255);
            }

            rgbaBuf[ot*(imageZ*imageY*imageX) + oz*(imageY*imageX) + oy*(imageX) + ox] = rgba;
        }
}
