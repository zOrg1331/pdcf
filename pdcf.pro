
#QT += gui
#CONFIG += gui

# to disable graphical user interface uncomment here
QT -= gui
CONFIG -= gui

TARGET = pdcf
CONFIG -= app_bundle
TEMPLATE = app
QMAKE_CXXFLAGS += -O3
SOURCES += main.cpp \
    pdcf.cpp \
    ls.cpp \
    common_math_tools.cpp \
    calccrow.cpp \
    calcrblock.cpp \
    pdcf_shell.cpp \
    pdcfcalc.cpp

gui {
    SOURCES += mainwindow.cpp
}

HEADERS += pdcf.h \
    ls.h \
    common_math_tools.h \
    calccrow.h \
    calcrblock.h \
    pdcf_shell.h \
    pdcfcalc.h

gui {
    HEADERS += mainwindow.h
}

gui {
    FORMS += mainwindow.ui
}

LIBS += -L./ -lboost_program_options
