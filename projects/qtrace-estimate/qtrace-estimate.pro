include(../base.pro)

TARGET = qtrace-estimate
CONFIG(debug, debug|release) {
	TARGET = $$join(TARGET,,,_debug)
}

RESOURCES += ../../src/qtrace.qrc

HEADERS += \
	../../src/Estimator.h \

SOURCES += \
	../../src/Estimator.cpp \
	../../src/qtrace-estimate.cpp \
