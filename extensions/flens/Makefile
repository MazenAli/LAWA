include $(LAWA_HOME)/Makefile.common

SOURCES  := $(wildcard *.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)

DYLIB    = libextensionsflens.$(DYLIB_EXT)
STATLIB  = libextensionsflens.$(STATLIB_EXT)

ifeq ($(TARGET), Darwin)
	LIB = $(DYLIB)
else
	LIB = $(STATLIB)
endif

all: $(LIB)

$(DYLIB): $(OBJECTS)
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

$(STATLIB): $(OBJECTS)
	ar rus $@ *.o
	$(CP) $(STATLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB) $(STATLIB)
	$(RM) -f *.tmp
	$(RMDIR) *.dSYM

# automatic dependencies -------------------------------------------------------
%.d: %.cpp
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
