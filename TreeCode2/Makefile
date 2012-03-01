include common.mk

DOC_CONF = doxy.conf
TARGET = tree

MAIN_SRCS = main.cpp $(SRCS)
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)


all: build

build: $(MAIN_OBJS)
	$(CC) $(LDFLAGS) -L. $(MAIN_OBJS) -o $(TARGET)

clean:
	rm -f $(OBJS)
	rm -f main.o
	rm -f ${TARGET}

.PHONY: doc
doc:
	doxygen $(DOC_CONF)



