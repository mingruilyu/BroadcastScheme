CXX = g++
CXXFLAGS = -g -W -Wall -pedantic -std=c++0x

SRCS = $(shell ls *.cpp)

APP = prog
OBJS = $(SRCS:.cpp=.o)

INCLUDES = -I../headers
all: $(APP)

$(APP) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP) $(OBJS)

%.o : %.cpp
	$(CXX) $< $(CXXFLAGS) $(INCLUDES) -c -o $@

clean:
	rm -rf *.o $(APP)
