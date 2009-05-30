TEMPLATE = app

CONFIG += debug_and_release console

DEPENDPATH += .
INCLUDEPATH += .

macx {
	CONFIG -= app_bundle
	CONFIG += x86
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

QT = core gui xml svg
