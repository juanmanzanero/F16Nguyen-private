$(shell mkdir -p ../lib/Release)


CXX=/opt_cluster/gcc/gcc-8.2.0/bin/g++ -std=c++2a
INCLUDE=-I../tr7/include -I../include -I/home/c01973/apps/include/lion/thirdparty/include
SRC=../src/F16_Nguyen_clib.cpp
OUT=../lib/Release/libF16_Nguyen_clib.so
FLAGS=-Ofast -flto -fno-fast-math -fmax-errors=10
DEFS=-DF16_NGUYEN_CLIB_EXPORTS -DTR7_WITH_AD


.PHONY:all
all:
	$(CXX) -fPIC -shared -fvisibility=hidden $(FLAGS) $(DEFS) $(INCLUDE) $(SRC) -o $(OUT)