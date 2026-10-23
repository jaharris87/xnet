# Compile external callers with the configuration exported by the production build.
# First run make -C ../.. BUILD_NAME=library-test libxnet net_setup.
XNET_CONFIG ?= ../../build/library-test/xnet-config.mk
include $(XNET_CONFIG)

BUILD_DIR ?= build
PROGRAMS := $(BUILD_DIR)/api_consumer $(BUILD_DIR)/legacy_reference

.PHONY: all test
all: $(PROGRAMS)

$(BUILD_DIR):
	mkdir -p '$@'

$(BUILD_DIR)/consumer_fixture.o: consumer_fixture.F90 $(XNET_CONFIG) | $(BUILD_DIR)
	$(XNET_FC) $(XNET_FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c '$<' -o '$@'

$(BUILD_DIR)/api_consumer.o: api_consumer.F90 $(BUILD_DIR)/consumer_fixture.o $(XNET_CONFIG)
	$(XNET_FC) $(XNET_FFLAGS) $(XNET_MODFLAGS) -I$(BUILD_DIR) -c '$<' -o '$@'

$(BUILD_DIR)/legacy_reference.o: legacy_reference.F90 $(BUILD_DIR)/consumer_fixture.o $(XNET_CONFIG)
	$(XNET_FC) $(XNET_FFLAGS) $(XNET_MODFLAGS) -I$(BUILD_DIR) -c '$<' -o '$@'

$(BUILD_DIR)/api_consumer: $(BUILD_DIR)/api_consumer.o $(BUILD_DIR)/consumer_fixture.o $(XNET_ARCHIVE)
	$(XNET_LDR) $(XNET_LDFLAGS) -o '$@' $(filter %.o,$^) $(XNET_ARCHIVE) $(XNET_LIBS)

$(BUILD_DIR)/legacy_reference: $(BUILD_DIR)/legacy_reference.o $(BUILD_DIR)/consumer_fixture.o $(XNET_ARCHIVE)
	$(XNET_LDR) $(XNET_LDFLAGS) -o '$@' $(filter %.o,$^) $(XNET_ARCHIVE) $(XNET_LIBS)

test: all
	python3 run_library_tests.py --api $(BUILD_DIR)/api_consumer --reference $(BUILD_DIR)/legacy_reference \
	  --net-setup $(dir $(XNET_ARCHIVE))../bin/net_setup --work-dir $(BUILD_DIR)/runs
