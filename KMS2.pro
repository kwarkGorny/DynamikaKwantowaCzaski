TEMPLATE = app
CONFIG += console c++1z
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O2 -DNDEBUG -fno-exceptions -flto -ffunction-sections -fdata-sections -finline-limit=150 -march=native -ffast-math

QMAKE_CXXFLAGS += -Winline -std=c++17

SOURCES += main.cpp

HEADERS += \
    WaveFunction.h
