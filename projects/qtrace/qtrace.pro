include(../base.pro)

TARGET = qtrace
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

RESOURCES += ../../src/qtrace.qrc

HEADERS += \
	../../src/Quantifier.h \

SOURCES += \
	../../src/Quantifier.cpp \
	../../src/qtrace.cpp \
