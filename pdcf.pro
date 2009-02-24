# -------------------------------------------------
# Project created by QtCreator 2009-02-03T16:35:17
# -------------------------------------------------
QT += gui
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
    mainwindow.cpp
HEADERS += pdcf.h \
    ls.h \
    common_math_tools.h \
    calccrow.h \
    calcrblock.h \
    mainwindow.h
FORMS += mainwindow.ui
