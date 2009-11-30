include(../base.pro)

TARGET = qtrace
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

RESOURCES += ../../src/qtrace.qrc

HEADERS += \
    ../../src/IsotopeEnvelope.h \
	../../src/Quantifier.h \
	../../src/Tango.h \

SOURCES += \
    ../../src/IsotopeEnvelope.cpp \
	../../src/Quantifier.cpp \
	../../src/qtrace.cpp \
