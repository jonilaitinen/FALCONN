CXX = g++
CXXFLAGS=-std=gnu++11 -DNDEBUG -Wall -Wextra -march=native -O3 -fPIC \
-I ../external/eigen \
-I ../src/include

LDFLAGS = -shared  # linking flags
RM = rm -f  # rm command
TARGET_LIB = libfalconn.so # target lib
TARGET_TEST = test

SRCS = FalconnSearchC.cpp FalconnSearch.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): $(OBJS)
	$(CXX) ${LDFLAGS} -o $@ $^

$(TARGET_TEST): $(OBJS)
	$(CXX) -o $@ $^

$(SRCS:.cpp=.d):%.d:%.cpp
	$(CXX) $(CXXFLAGS) -MM $< >$@

include $(SRCS:.cpp=.d)

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB} ${OBJS} $(SRCS:.cpp=.d)