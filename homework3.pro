TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    mat4.cpp \
    tiny_obj_loader.cc \
    vec4.cpp \
    others.cpp

HEADERS += \
    mat4.h \
    tiny_obj_loader.h \
    vec4.h \
    others.h

RESOURCES +=
