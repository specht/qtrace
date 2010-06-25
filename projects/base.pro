TEMPLATE = app

CONFIG += debug_and_release console

DEPENDPATH += .
INCLUDEPATH += .

macx {
    LIBPATH += /Users/michael/programming/ext/lib
    INCLUDEPATH += /Users/michael/programming/ext/include
}

macx {
    CONFIG -= app_bundle
    CONFIG += x86 ppc
    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.4
}

unix {
    CONFIG += static
    QMAKE_LFLAGS += -L/usr/local/lib
    LIBS += -ldl -lrt
    message("STATIC LINKAGE")
}

LIBS += -lptb -lz -lbz2 -lquazip

CONFIG(debug, debug|release) {
    OBJECTS_DIR = ../../obj/debug/
    MOC_DIR = ../../obj/debug/
    RCC_DIR = ../../obj/debug/
}
else {
    OBJECTS_DIR = ../../obj/release/
    MOC_DIR = ../../obj/release/
    RCC_DIR = ../../obj/release/
}

DESTDIR = ../../

QT = core xml

HEADERS += \
    ../../src/QuantifierBase.h \
    ../../src/Tango.h \

SOURCES += \
    ../../src/QuantifierBase.cpp \
    ../../src/Tango.cpp \
