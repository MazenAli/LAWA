include $(LAWA_HOME)/Makefile.common

DYLIB    = libextensionssparsegrid.$(DYLIB_EXT)

SOURCES  := sandia_rules.cpp sparse_grid_mixed.cpp sparse_grid_open.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)

DYLIB    = libextensionssparsegrid.$(DYLIB_EXT)
STATLIB  = libextensionssparsegrid.$(STATLIB_EXT)

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
