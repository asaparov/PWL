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
EXTRACT_WIKT_MORPHOLOGY_EN_CPP_SRCS=extract_wikt_morphology_en.cpp
EXTRACT_WIKT_MORPHOLOGY_EN_DBG_OBJS=$(EXTRACT_WIKT_MORPHOLOGY_EN_CPP_SRCS:.cpp=.debug.o)
EXTRACT_WIKT_MORPHOLOGY_EN_OBJS=$(EXTRACT_WIKT_MORPHOLOGY_EN_CPP_SRCS:.cpp=.release.o)
SET_REASONING_TEST_CPP_SRCS=set_reasoning_test.cpp
SET_REASONING_TEST_DBG_OBJS=$(SET_REASONING_TEST_CPP_SRCS:.cpp=.debug.o)
SET_REASONING_TEST_OBJS=$(SET_REASONING_TEST_CPP_SRCS:.cpp=.release.o)


#
# Compile and link options
#

CPP=g++
cc-option = $(shell $(CPP) -Werror $(1) -c -x c /dev/null -o /dev/null 2>/dev/null; echo $$?)

LIBRARY_PKG_LIBS=
PKG_LIBS=-pthread -lssl -lcrypto
NO_AS_NEEDED=-Wl,--no-as-needed
ifeq ($(call cc-option, $(NO_AS_NEEDED)),0)
	PKG_LIBS += $(NO_AS_NEEDED)
endif
GLIBC := $(word 2,$(shell getconf GNU_LIBC_VERSION 2>/dev/null))
ifeq "$(.SHELLSTATUS)" "0"
	GLIBC_HAS_RT := $(shell expr $(GLIBC) \>= 2.17)
	ifeq "$(GLIBC_HAS_RT)" "0"
		LIBRARY_PKG_LIBS += -lrt
		PKG_LIBS += -lrt
	endif
endif

WARNING_FLAGS=-Wall -Wpedantic
override CPPFLAGS_DBG += $(WARNING_FLAGS) -I. -g -march=native -mtune=native -std=c++11
override CPPFLAGS += $(WARNING_FLAGS) -I. -Ofast -fno-finite-math-only -DNDEBUG -march=native -mtune=native -std=c++11 -fno-stack-protector
override LDFLAGS_DBG += -g $(LIB_PATHS) $(PKG_LIBS)
override LDFLAGS += $(LIB_PATHS) -fwhole-program $(PKG_LIBS)


#
# GNU Make: targets that don't build files
#

.PHONY: all debug clean distclean

#
# Make targets
#

all: executive_test extract_wikt_morphology_en set_reasoning_test

debug: executive_test_dbg extract_wikt_morphology_en_dbg set_reasoning_test_dbg

-include $(EXECUTIVE_TEST_OBJS:.release.o=.release.d)
-include $(EXECUTIVE_TEST_DBG_OBJS:.debug.o=.debug.d)
-include $(EXTRACT_WIKT_MORPHOLOGY_EN_OBJS:.release.o=.release.d)
-include $(EXTRACT_WIKT_MORPHOLOGY_EN_DBG_OBJS:.debug.o=.debug.d)
-include $(SET_REASONING_TEST_OBJS:.release.o=.release.d)
-include $(SET_REASONING_TEST_DBG_OBJS:.debug.o=.debug.d)

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

executive_test: $(LIBS) $(EXECUTIVE_TEST_OBJS)
		$(CPP) -o executive_test $(CPPFLAGS) $(EXECUTIVE_TEST_OBJS) $(LDFLAGS)

executive_test_dbg: $(LIBS) $(EXECUTIVE_TEST_DBG_OBJS)
		$(CPP) -o executive_test_dbg $(CPPFLAGS_DBG) $(EXECUTIVE_TEST_DBG_OBJS) $(LDFLAGS_DBG)

extract_wikt_morphology_en: $(LIBS) $(EXTRACT_WIKT_MORPHOLOGY_EN_OBJS)
		$(CPP) -o extract_wikt_morphology_en $(CPPFLAGS) $(EXTRACT_WIKT_MORPHOLOGY_EN_OBJS) $(LDFLAGS)

extract_wikt_morphology_en_dbg: $(LIBS) $(EXTRACT_WIKT_MORPHOLOGY_EN_DBG_OBJS)
		$(CPP) -o extract_wikt_morphology_en_dbg $(CPPFLAGS_DBG) $(EXTRACT_WIKT_MORPHOLOGY_EN_DBG_OBJS) $(LDFLAGS_DBG)

set_reasoning_test: $(LIBS) $(SET_REASONING_TEST_OBJS)
		$(CPP) -o set_reasoning_test $(CPPFLAGS) $(SET_REASONING_TEST_OBJS) $(LDFLAGS)

set_reasoning_test_dbg: $(LIBS) $(SET_REASONING_TEST_DBG_OBJS)
		$(CPP) -o set_reasoning_test_dbg $(CPPFLAGS_DBG) $(SET_REASONING_TEST_DBG_OBJS) $(LDFLAGS_DBG)

clean:
	    ${RM} -f *.o */*.o */*/*.o *.d */*.d */*/*.d executive_test executive_test.exe executive_test_dbg executive_test_dbg.exe extract_wikt_morphology_en extract_wikt_morphology_en_dbg extract_wikt_morphology_en.exe extract_wikt_morphology_en_dbg.exe set_reasoning_test set_reasoning_test_dbg set_reasoning_test.exe set_reasoning_test_dbg.exe $(LIBS)
