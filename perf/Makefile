SINCWIN=sincwin_ref sincwin_linear sincwin_linear_polyphase sincwin_quadratic sincwin_quadratic_polyphase

SOURCES=$(SINCWIN:%=%.c) perf.cpp
EXEC=$(SINCWIN:%=perf_%)

FLAGS=-Wall -O3  -ffast-math    # -march=native
CFLAGS=
CXXFLAGS=-std=c++11
CPPFLAGS=$(CFLAGS) $(FLAGS)

OBJECTS=$(patsubst %.c,$(PREFIX)%$(SUFFIX).o,\
        $(patsubst %.cpp,$(PREFIX)%$(SUFFIX).o,\
$(SOURCES)))
DEPENDS=$(patsubst %.c,$(PREFIX)%$(SUFFIX).d,\
        $(patsubst %.cpp,$(PREFIX)%$(SUFFIX).d,\
        $(filter-out %.o,""\
$(SOURCES))))

all: $(EXEC)
ifneq "$(MAKECMDGOALS)" "clean"
  -include $(DEPENDS)
endif

clean:
	rm -f $(EXEC) $(OBJECTS) $(DEPENDS)

%.d: %.c
	@$(CC) $(CPPFLAGS) -MM -MT"$@" -MT"$*.o" -o $@ $<  2> /dev/null

%.d: %.cpp
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MT"$@" -MT"$*.o" -o $@ $<  2> /dev/null

$(EXEC): perf_%: %.o perf.o
	$(CXX) -o $@ $^ $(LDFLAGS)

perf: $(filter-out perf_sincwin_ref,$(EXEC))
	@for i in $^; do echo $$i; ./$$i; done

