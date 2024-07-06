SOURCES=resampler.cpp
EXEC=test1 cli1

CFLAGS=-Wall -O3 -ffast-math
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
	rm -f $(EXEC) $(EXEC:%=%.o) $(OBJECTS) $(DEPENDS) spectrogram_sweep.png

%.d: %.c
	@$(CC) $(CPPFLAGS) -MM -MT"$@" -MT"$*.o" -o $@ $<  2> /dev/null

%.d: %.cpp
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MT"$@" -MT"$*.o" -o $@ $<  2> /dev/null

$(EXEC): %: %.o $(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: spectrograms
spectrograms: spectrogram_sweep.png

spectrogram_sweep.png: cli1
	sox -r 96000 -n -t f32 - synth 8 sine 0+48000 gain -6.0206 | ./cli1 44100/96000 | sox -t f32 -r 44100 -c 1 - -n spectrogram -c "" -w Kaiser -o $@

