include(../base.pro)

TARGET = qtrace-estimate
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

RESOURCES += ../../src/qtrace.qrc

HEADERS += \
	../../src/Quantifier.h \
	../../src/Tango.h \

SOURCES += \
	../../src/Quantifier.cpp \
	../../src/qtrace-estimate.cpp \
