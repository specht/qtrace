TEMPLATE = app

CONFIG += debug_and_release console

DEPENDPATH += .
INCLUDEPATH += .
LIBS += -lbz2

macx {
	CONFIG -= app_bundle
	CONFIG += x86
}

win32 {
	LIBS += -lzdll
} else {
	LIBS += -lz
}

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

QT = core gui xml svg

HEADERS += \
	MzDataHandler.h \
	MzMlHandler.h \
	MzXmlHandler.h \
	RefPtr.h \
	ScanIterator.h \
	XmlHandler.h \
	ZipFileOrNot.h \

SOURCES += \
	MzDataHandler.cpp \
	MzMlHandler.cpp \
	MzXmlHandler.cpp \
	ScanIterator.cpp \
	XmlHandler.cpp \
	ZipFileOrNot.cpp \

# QuaZip

DEPENDPATH += ../../src/
INCLUDEPATH += ../../src/quazip-0.2.2/quazip/

# Input
HEADERS += quazip-0.2.2/quazip/crypt.h \
           quazip-0.2.2/quazip/ioapi.h \
           quazip-0.2.2/quazip/quazip.h \
           quazip-0.2.2/quazip/quazipfile.h \
           quazip-0.2.2/quazip/quazipfileinfo.h \
           quazip-0.2.2/quazip/quazipnewinfo.h \
           quazip-0.2.2/quazip/unzip.h \
           quazip-0.2.2/quazip/zip.h

SOURCES += quazip-0.2.2/quazip/ioapi.c \
           quazip-0.2.2/quazip/quazip.cpp \
           quazip-0.2.2/quazip/quazipfile.cpp \
           quazip-0.2.2/quazip/quazipnewinfo.cpp \
           quazip-0.2.2/quazip/unzip.c \
           quazip-0.2.2/quazip/zip.c
           
# KFilterBase

INCLUDEPATH += ../../src/kfilterbase/

# Input

HEADERS += \
	kfilterbase/kbzip2filter.h \
	kfilterbase/kde_config.h \
	kfilterbase/kfilterbase.h \
	kfilterbase/kfilterdev.h \
	kfilterbase/kgzipfilter.h \

SOURCES += \
	kfilterbase/kbzip2filter.cpp \
	kfilterbase/kfilterbase.cpp \
	kfilterbase/kfilterdev.cpp \
	kfilterbase/kgzipfilter.cpp	 \
