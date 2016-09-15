TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    oppg_c.cpp
LIBS += -llapack -lblas -larmadillo
