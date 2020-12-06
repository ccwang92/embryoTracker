QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0
DEFINES += USE_Qt5 # for vaa3d functions

SOURCES += emt_glwidget.cpp \
    data_importer.cpp \
    main.cpp \
    mainwindow.cpp \
    my4dimage.cpp \
    src_3rd/basic_c_fun/basic_4dimage.cpp \
    src_3rd/basic_c_fun/basic_4dimage_create.cpp \
    src_3rd/basic_c_fun/imageio_mylib.cpp \
    src_3rd/basic_c_fun/mg_image_lib.cpp \
    src_3rd/basic_c_fun/mg_utilities.cpp \
    src_3rd/basic_c_fun/stackutil.cpp \
    src_3rd/basic_c_fun/v3d_message.cpp \
    src_3rd/io/io_bioformats.cpp

HEADERS += emt_glwidget.h \
    data_importer.h \
    mainwindow.h \
    my4dimage.h \
    src_3rd/basic_c_fun/basic_4dimage.h \
    src_3rd/basic_c_fun/imageio_mylib.h \
    src_3rd/basic_c_fun/mg_image_lib.h \
    src_3rd/basic_c_fun/mg_utilities.h \
    src_3rd/basic_c_fun/stackutil.h \
    src_3rd/basic_c_fun/v3d_basicdatatype.h \
    src_3rd/basic_c_fun/v3d_message.h \
    src_3rd/io/io_bioformats.h \
    src_3rd/v3d/import_images_tool_dialog.h

#FORMS += v3d/import_images_tool.ui

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
    src_3rd/v3d/import_images_tool.ui



LIBS += -L$$PWD/src_3rd/common_lib/lib
LIBS += -lm -lv3dtiff # for tiff read
#LIBS += -L$$PWD/src_3rd/jba/c++
#LIBS += -lv3dnewmat  # for new image create
LIBS += -L$$PWD/src_3rd/common_lib/src_packages/mylib_tiff/ -lmylib # for usage of self-deisgned tiff functions
#LIBS += -L$$PWD/src_3rd/common_lib/lib_unix64  -lteem  -lbz2 -lz  -lGLU #for nrrd support

#INCLUDEPATH += $$PWD/src_3rd/common_lib/include
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
