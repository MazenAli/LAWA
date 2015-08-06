include $(LAWA_HOME)/Makefile.common

SUBDIRS := math aux methods
SOURCES  := $(wildcard *.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(SOURCES:.cpp=.d)

DYLIB    = liblawa.$(DYLIB_EXT)
STATLIB  = liblawa.$(STATLIB_EXT)

ifeq ($(TARGET), Darwin)
	LIB = $(DYLIB)
else
	LIB = $(STATLIB)
endif

all: $(LIB)

$(DYLIB): $(OBJECTS)
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
	$(CXX) $(CXXDYLIB) $(CXXFLAGS) $(LDFLAGS) -o $(DYLIB) *.o aux/*.o math/*.o methods/adaptive/datastructures/*.o $(LDLIBS)
	$(CP) $(DYLIB) $(LAWA_HOME)

$(STATLIB): $(OBJECTS)
	$(SILENT) for i in $(SUBDIRS); \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i all; \
	done;
	ar rus $@ *.o aux/*.o math/*.o methods/adaptive/datastructures/*.o
	$(CP) $(STATLIB) $(LAWA_HOME)

clean:
	$(RM) *.o *.d $(DYLIB) $(STATLIB)
	$(SILENT) for i in $(SUBDIRS) ; \
		do echo ""; echo "processing dir $$i"; \
		$(MAKE) -C $$i clean; \
	done;
	$(RM) -f *.tmp

# automatic dependencies -------------------------------------------------------
%.d: %.cpp
	$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< \
		| sed 's/$*.o/& $(subst /,\/,$@)/g' > $@

-include $(DEPS)
