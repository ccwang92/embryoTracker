QT       += core opengl #gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
#DEFINES += USE_Qt5 # for vaa3d functions

SOURCES += data_importer.cpp \
    cellsegmentation/cellsegment_main.cpp \
    cellsegmentation/img_basic_proc.cpp \
    cellsegmentation/maxflow_bk/graph.cpp \
    cellsegmentation/maxflow_bk/maxflow.cpp \
    cellsegmentation/synquant_simple.cpp \
    celltracking/CINDA/src_c/cinda_funcs.c \
    celltracking/celltracking_main.cpp \
    main.cpp \
    mainwindow.cpp \
    my4dimage.cpp \
    myglwidget.cpp \
#    emt_glwidget.cpp \
#    renderer_emt.cpp \
#    src_3rd/3drenderer/GLee2glew.c \
#    src_3rd/3drenderer/renderer.cpp \
#    src_3rd/3drenderer/renderer_gl2.cpp \
#    src_3rd/3drenderer/glsl_r.cpp \
    raycasting/mesh.cpp \
    raycasting/raycastcanvas.cpp \
    raycasting/raycastvolume.cpp \
    raycasting/trackball.cpp \
    src_3rd/basic_c_fun/basic_4dimage.cpp \
    src_3rd/basic_c_fun/basic_4dimage_create.cpp \
    src_3rd/basic_c_fun/imageio_mylib.cpp \
    src_3rd/basic_c_fun/mg_image_lib.cpp \
    src_3rd/basic_c_fun/mg_utilities.cpp \
    src_3rd/basic_c_fun/stackutil.cpp \
    src_3rd/basic_c_fun/v3d_message.cpp \
    src_3rd/cellseg/FL_morphology.cpp \
    src_3rd/io/io_bioformats.cpp

HEADERS += data_importer.h \
    cellsegmentation/cc3d.hpp \
    cellsegmentation/cellsegment_main.h \
    cellsegmentation/img_basic_proc.h \
    cellsegmentation/img_basic_proc_declare.h \
    cellsegmentation/maxflow_bk/block.h \
    cellsegmentation/maxflow_bk/graph.h \
    cellsegmentation/synquant_simple.h \
    cellsegmentation/types_define.h \
    cellsegmentation/vol_basic_proc.hpp \
    celltracking/celltracking_main.h \
    celltracking/dt3d.hpp \
    dialog_import_im_sequence.h \
    mainwindow.h \
    my4dimage.h \
    myglwidget.h \
#    emt_glwidget.h \
#    renderer_emt.h \
#    src_3rd/3drenderer/GLee2glew.h \
#    src_3rd/3drenderer/renderer.h \
#    src_3rd/3drenderer/renderer_gl1.h \
#    src_3rd/3drenderer/renderer_gl2.h \
#    src_3rd/3drenderer/glsl_r.h \
    raycasting/mesh.h \
    raycasting/raycastcanvas.h \
    raycasting/raycastvolume.h \
    raycasting/trackball.h \
    src_3rd/basic_c_fun/basic_4dimage.h \
    src_3rd/basic_c_fun/color_xyz.h \
    src_3rd/basic_c_fun/imageio_mylib.h \
    src_3rd/basic_c_fun/mg_image_lib.h \
    src_3rd/basic_c_fun/mg_utilities.h \
    src_3rd/basic_c_fun/stackutil.h \
    src_3rd/basic_c_fun/v3d_basicdatatype.h \
    src_3rd/basic_c_fun/v3d_message.h \
    src_3rd/basic_c_fun/volimg_proc.h \
    src_3rd/basic_c_fun/volimg_proc_declare.h \
    src_3rd/cellseg/FL_accessory.h \
    src_3rd/cellseg/FL_adaptiveThreshold3D.h \
    src_3rd/cellseg/FL_bwdist.h \
    src_3rd/cellseg/FL_bwlabel2D3D.h \
    src_3rd/cellseg/FL_cellSegmentation3D.h \
    src_3rd/cellseg/FL_coordDefinition.h \
    src_3rd/cellseg/FL_defType.h \
    src_3rd/cellseg/FL_distanceTransform3D.h \
    src_3rd/cellseg/FL_downSample3D.h \
    src_3rd/cellseg/FL_filter3D.h \
    src_3rd/cellseg/FL_gvfCellSeg.h \
    src_3rd/cellseg/FL_interpolateCoordCubic3D.h \
    src_3rd/cellseg/FL_interpolateCoordLinear3D.h \
    src_3rd/cellseg/FL_main_brainseg.h \
    src_3rd/cellseg/FL_morphology.h \
    src_3rd/cellseg/FL_neighborhood.h \
    src_3rd/cellseg/FL_neighborhoodWalker.h \
    src_3rd/cellseg/FL_queue.h \
    src_3rd/cellseg/FL_regionProps.h \
    src_3rd/cellseg/FL_sort.h \
    src_3rd/cellseg/FL_sort2.h \
    src_3rd/cellseg/FL_threshold.h \
    src_3rd/cellseg/FL_unionFind.h \
    src_3rd/cellseg/FL_upSample3D.h \
    src_3rd/cellseg/FL_volimgProcLib.h \
    src_3rd/cellseg/FL_watershed_vs.h \
    src_3rd/cellseg/nrutil.h \
    src_3rd/io/io_bioformats.h \
#    src_3rd/v3d/import_images_tool_dialog.h

#FORMS += v3d/import_images_tool.ui


## adding include paths and libs related to opencv
#include(opencv.qrc)
INCLUDEPATH += /usr/local/include/opencv4
LIBS += -L/usr/local/lib -lopencv_imgproc -lopencv_core -lopencv_highgui -lopencv_imgcodecs

# include(boost.qrc)
#INCLUDEPATH += /usr/include/
#LIBS += -L/usr/include/boost/bin.v2/libs/ -lboost_system -lboost_filesystem -lboost_asio
#include(VTK.qrc)
#INCLUDEPATH += /usr/local/include/ITK-5.1
#DEPENDPATH += /usr/local/include/ITK-5.1
#LIBS += -L/usr/local/lib -lITKBiasCorrection-5.1 -lITKBioCell-5.1 -lITKCommon-5.1 -lITKDICOMParser-5.1 -litkdouble-conversion-5.1 -lITKEXPAT-5.1 -lITKFEM-5.1 -litkgdcmCommon-5.1 -litkgdcmDICT-5.1 -litkgdcmDSED-5.1\
#        -litkgdcmIOD-5.1 -litkgdcmjpeg8-5.1 -litkgdcmjpeg12-5.1 -litkgdcmjpeg16-5.1 -litkgdcmMSFF-5.1 -lITKgiftiio-5.1 -litkhdf5_cpp-5.1 -litkhdf5-5.1 -lITKIOBioRad-5.1 -lITKIOBMP-5.1\
#        -lITKIOCSV-5.1 -lITKIOGDCM-5.1 -lITKIOGE-5.1 -lITKIOGIPL-5.1 -lITKIOHDF5-5.1 -lITKIOImageBase-5.1 -lITKIOIPL-5.1 -lITKIOJPEG-5.1 -lITKIOLSM-5.1 -lITKIOMesh-5.1 -lITKIOMeta-5.1 -lITKIOMRC-5.1\
#        -lITKIONIFTI-5.1 -lITKIONRRD-5.1 -lITKIOPNG-5.1 -lITKIOSiemens-5.1 -lITKIOSpatialObjects-5.1 -lITKIOStimulate-5.1 -lITKIOTIFF-5.1 -lITKIOTransformBase-5.1 -lITKIOTransformHDF5-5.1 -lITKIOTransformInsightLegacy-5.1\
#        -lITKIOTransformMatlab-5.1 -lITKIOVTK-5.1 -lITKIOXML-5.1 -litkjpeg-5.1 -lITKKLMRegionGrowing-5.1 -lITKLabelMap-5.1 -lITKMesh-5.1 -lITKMetaIO-5.1 -litkNetlibSlatec-5.1 -lITKniftiio-5.1\
#        -lITKNrrdIO-5.1 -litkopenjpeg-5.1 -lITKOptimizers-5.1 -lITKPath-5.1 -litkpng-5.1 -lITKPolynomials-5.1 -lITKQuadEdgeMesh-5.1 -lITKSpatialObjects-5.1 -lITKStatistics-5.1\
#        -litktiff-5.1 -litkv3p_lsqr-5.1 -litkv3p_netlib-5.1 -litkvcl-5.1 -lITKVideoCore-5.1 -lITKVideoIO-5.1 -litkvnl_algo-5.1 -litkvnl-5.1 -lITKVNLInstantiation-5.1 -lITKVTK-5.1\
#        -lITKVtkGlue-5.1 -lITKWatersheds-5.1 -litkzlib-5.1 -lITKznz-5.1 -lITKVideoBridgeOpenCV-5.1
#LIBS += src_3d/common_lib/src_packages/mylib_tiff -lmylib

#INCLUDEPATH += src_3rd/common_lib/include/glew
#INCLUDEPATH += src_3rd/v3d
#INCLUDEPATH += src_3rd/3drenderer
#INCLUDEPATH += src_3rd/common_lib/include
#INCLUDEPATH += src_3rd/basic_c_fun
#INCLUDEPATH = $$unique(INCLUDEPATH)
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

FORMS += \
    dialog_import_im_sequence.ui #\
    mainwindow.ui
#    src_3rd/v3d/import_images_tool.ui



LIBS += -L$$PWD/src_3rd/common_lib/lib
LIBS += -lm -lv3dtiff # for tiff read
#LIBS += -L$$PWD/src_3rd/jba/c++
#LIBS += -lv3dnewmat  # for new image create
LIBS += -L$$PWD/src_3rd/common_lib/src_packages/mylib_tiff/ -lmylib # for usage of self-deisgned tiff functions
#LIBS += -L$$PWD/src_3rd/common_lib/lib_unix64  -lteem  -lbz2 -lz  -lGLU #for nrrd support
#LIBS += -lglut -lGLU # for GL support

#INCLUDEPATH += $$PWD/src_3rd/common_lib/include # for glew folder
#INCLUDEPATH += $$PWD/src_3rd/common_lib/include/glew/
#INCLUDEPATH += $$PWD/src_3rd/basic_c_fun


#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/src_3rd/common_lib/src_packages/mylib_tiff/release/ -lmylib
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/src_3rd/common_lib/src_packages/mylib_tiff/debug/ -lmylib
#else:unix: LIBS += -L$$PWD/src_3rd/common_lib/src_packages/mylib_tiff/ -lmylib



#INCLUDEPATH += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff
#DEPENDPATH += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff/release/libmylib.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff/debug/libmylib.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff/release/mylib.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff/debug/mylib.lib
#else:unix: PRE_TARGETDEPS += $$PWD/src_3rd/common_lib/src_packages/mylib_tiff/libmylib.a

RESOURCES += \
    qdarkstyle/style.qrc \
    raycasting/resources.qrc
