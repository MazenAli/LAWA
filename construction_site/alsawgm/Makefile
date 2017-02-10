include $(LAWA_HOME)/Makefile.common

SOURCEFILES  = $(wildcard *.cpp)
TARGETS      = $(patsubst %.cpp,%,$(SOURCEFILES))

LDFLAGS += -L$(LAWA_HOME)
LDLIBS  := -llawa -lextensionsflens -lextensionssparsegrid $(LDLIBS)

ifeq ($(TARGET), Darwin)
	EX_LDFLAGS += -bind_at_load
endif

ifdef GMP_HOME
	  CXXFLAGS += -DMPFR -I$(MPFR_HOME)/include
	  LDFLAGS += -L$(MPFR_HOME)/lib
	  LDLIBS  += -lgmpxx -lgmp
endif

all: $(TARGETS)

clean:
		$(RM) *.o
			$(RM) $(TARGETS)
				$(RMDIR) *.dSYM
