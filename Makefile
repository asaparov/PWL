#
# Unix/Linux makefile
# Abulhair Saparov
#

#
# List of source files
#

EXECUTIVE_TEST_CPP_SRCS=executive_test.cpp
EXECUTIVE_TEST_DBG_OBJS=$(EXECUTIVE_TEST_CPP_SRCS:.cpp=.debug.o)
EXECUTIVE_TEST_OBJS=$(EXECUTIVE_TEST_CPP_SRCS:.cpp=.release.o)


#
# Compile and link options
#

LIBRARY_PKG_LIBS=
PKG_LIBS=-Wl,--no-as-needed -lpthread
GLIBC := $(word 2,$(shell getconf GNU_LIBC_VERSION))
GLIBC_HAS_RT := $(shell expr $(GLIBC) \>= 2.17)
ifeq "$(GLIBC_HAS_RT)" "0"
	LIBRARY_PKG_LIBS += -lrt
	PKG_LIBS += -lrt
endif

CPP=g++
WARNING_FLAGS=-Wall -Wpedantic
override CPPFLAGS_DBG += $(WARNING_FLAGS) -I. -g -march=native -mtune=native -std=c++11 $(PKG_LIBS)
override CPPFLAGS += $(WARNING_FLAGS) -I. -Ofast -DNDEBUG -march=native -mtune=native -std=c++11 $(PKG_LIBS)
override LDFLAGS_DBG += -g $(LIB_PATHS)
override LDFLAGS += $(LIB_PATHS) -fwhole-program


#
# Compile command
#

-include $(EXECUTIVE_TEST_OBJS:.release.o=.release.d)
-include $(EXECUTIVE_TEST_DBG_OBJS:.debug.o=.debug.d)

define make_dependencies
	$(1) $(2) -c $(3).$(4) -o $(3).$(5).o
	$(1) -MM $(2) $(3).$(4) > $(3).$(5).d
	@mv -f $(3).$(5).d $(3).$(5).d.tmp
	@sed -e 's|.*:|$(3).$(5).o:|' < $(3).$(5).d.tmp > $(3).$(5).d
	@sed -e 's/.*://' -e 's/\\$$//' < $(3).$(5).d.tmp | fmt -1 | \
		sed -e 's/^ *//' -e 's/$$/:/' >> $(3).$(5).d
	@rm -f $(3).$(5).d.tmp
endef

%.release.o: %.cpp
	$(call make_dependencies,$(CPP),$(CPPFLAGS),$*,cpp,release)
%.release.pic.o: %.cpp
	$(call make_dependencies,$(CPP),$(CPPFLAGS),$*,cpp,release.pic)
%.debug.o: %.cpp
	$(call make_dependencies,$(CPP),$(CPPFLAGS_DBG),$*,cpp,debug)
%.debug.pic.o: %.cpp
	$(call make_dependencies,$(CPP),$(CPPFLAGS_DBG),$*,cpp,debug.pic)


#
# GNU Make: targets that don't build files
#

.PHONY: all debug clean distclean

#
# Make targets
#

all: executive_test

debug: executive_test_dbg

executive_test: $(LIBS) $(EXECUTIVE_TEST_OBJS)
		$(CPP) -o executive_test $(CPPFLAGS) $(LDFLAGS) $(EXECUTIVE_TEST_OBJS)

executive_test_dbg: $(LIBS) $(EXECUTIVE_TEST_DBG_OBJS)
		$(CPP) -o executive_test_dbg $(CPPFLAGS_DBG) $(LDFLAGS_DBG) $(EXECUTIVE_TEST_DBG_OBJS)

clean:
	    ${RM} -f *.o */*.o */*/*.o *.d */*.d */*/*.d executive_test executive_test.exe executive_test_dbg executive_test_dbg.exe $(LIBS)

distclean:  clean
	    ${RM} -f *~
