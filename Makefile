CXX ?= g++
CXXFLAGS ?= -O2 -g
BUILD ?= build
PROJECT_FLAGS = -std=c++20 -Wall -Wextra -Wpedantic -DCSV_ENABLE_THREADS=0
INCLUDES = -Iinclude -isystem third_party/csv-parser -isystem third_party/doctest

.PHONY: all test check verify-vendor sanitize
all: $(BUILD)/rnaseq

$(BUILD):
	mkdir -p $@

$(BUILD)/table.o: src/table.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/file_hash.o: src/file_hash.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/manifests.o: src/manifests.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/counts.o: src/counts.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/process.o: src/process.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/alignment.o: src/alignment.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_process.o: tests/unit/test_process.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/main.o: src/main.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_table.o: tests/unit/test_table.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_file_hash.o: tests/unit/test_file_hash.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_manifests.o: tests/unit/test_manifests.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_counts.o: tests/unit/test_counts.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/rnaseq: $(BUILD)/workflow.o $(BUILD)/process.o $(BUILD)/alignment.o $(BUILD)/main.o $(BUILD)/table.o $(BUILD)/file_hash.o $(BUILD)/manifests.o $(BUILD)/counts.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD)/test_table: $(BUILD)/process.o $(BUILD)/test_process.o $(BUILD)/test_table.o $(BUILD)/test_file_hash.o $(BUILD)/test_manifests.o $(BUILD)/test_counts.o $(BUILD)/table.o $(BUILD)/file_hash.o $(BUILD)/manifests.o $(BUILD)/counts.o
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
	$(MAKE) BUILD=build/sanitize CXXFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' LDFLAGS='-fsanitize=address,undefined' check

-include $(wildcard $(BUILD)/*.d)

$(BUILD)/workflow.o: src/workflow.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

# Export instead of interpolating user paths into recipe shell code.
export WORKFLOW_CONFIG
.PHONY: workflow
workflow: $(BUILD)/rnaseq
	@test -n "$$WORKFLOW_CONFIG" || { echo "Set WORKFLOW_CONFIG to a workflow TSV" >&2; exit 2; }
	$(BUILD)/rnaseq workflow-local "$$WORKFLOW_CONFIG"
