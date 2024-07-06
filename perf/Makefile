SINCWIN=sincwin_ref sincwin_nearest sincwin_nearest_polyphase sincwin_linear sincwin_linear_polyphase sincwin_quadratic sincwin_quadratic_polyphase

SOURCES=$(SINCWIN:%=%.c) perf.cpp approx.cpp
EXEC_PERF=$(SINCWIN:%=perf_%)
EXEC_APPROX=$(SINCWIN:%=approx_%)
EXEC=$(EXEC_PERF) $(EXEC_APPROX)

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

$(EXEC_PERF): perf_%: %.o perf.o
	$(CXX) -o $@ $^ $(LDFLAGS)

$(EXEC_APPROX): approx_%: %.o approx.o
	$(CXX) -o $@ $^ $(LDFLAGS)

perf: $(filter-out perf_sincwin_ref,$(EXEC_PERF))
	@for i in $^; do echo $$i; ./$$i; done

approx: $(EXEC_APPROX)
	@for i in $^; do echo $$i; ./$$i; done

