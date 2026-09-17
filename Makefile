CXX ?= g++
CXXFLAGS ?= -O2 -g
BUILD ?= build
PROJECT_FLAGS := -std=c++20 -Wall -Wextra -Wpedantic -DCSV_ENABLE_THREADS=0
INCLUDES := -Iinclude -isystem third_party/csv-parser -isystem third_party/doctest

CORE_SOURCES := src/table.cpp src/file_hash.cpp src/manifests.cpp src/counts.cpp src/process.cpp
CLI_SOURCES := src/main.cpp src/alignment.cpp src/workflow.cpp \
               src/workflow/config.cpp src/workflow/state.cpp src/workflow/execution.cpp
TEST_SOURCES := tests/unit/test_table.cpp tests/unit/test_file_hash.cpp \
                tests/unit/test_manifests.cpp tests/unit/test_counts.cpp tests/unit/test_process.cpp
CORE_OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(CORE_SOURCES))
CLI_OBJECTS := $(patsubst src/%.cpp,$(BUILD)/%.o,$(CLI_SOURCES))
TEST_OBJECTS := $(patsubst tests/unit/%.cpp,$(BUILD)/tests/%.o,$(TEST_SOURCES))
OBJECTS := $(CORE_OBJECTS) $(CLI_OBJECTS) $(TEST_OBJECTS)

.PHONY: all test check verify-vendor sanitize workflow
all: $(BUILD)/rnaseq

$(BUILD)/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/tests/%.o: tests/unit/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/rnaseq: $(CORE_OBJECTS) $(CLI_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD)/test_table: $(CORE_OBJECTS) $(TEST_OBJECTS)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

verify-vendor:
	sha256sum -c third_party/SHA256SUMS

test: $(BUILD)/test_table
	$(BUILD)/test_table

check: verify-vendor test $(BUILD)/rnaseq
	bash tests/integration/check_native.sh $(BUILD)/rnaseq
	bash tests/integration/check_counts.sh $(BUILD)/rnaseq
	bash tests/integration/check_p4_native.sh $(BUILD)/rnaseq
	bash tests/integration/check_p5_workflow.sh $(BUILD)/rnaseq

sanitize:
	$(MAKE) BUILD=$(BUILD)/sanitize CXXFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' LDFLAGS='-fsanitize=address,undefined' check

# Export instead of interpolating user paths into recipe shell code.
export WORKFLOW_CONFIG
workflow: $(BUILD)/rnaseq
	@test -n "$$WORKFLOW_CONFIG" || { echo "Set WORKFLOW_CONFIG to a workflow TSV" >&2; exit 2; }
	$(BUILD)/rnaseq workflow-local "$$WORKFLOW_CONFIG"

-include $(OBJECTS:.o=.d)
