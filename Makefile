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

$(BUILD)/main.o: src/main.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/test_table.o: tests/unit/test_table.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PROJECT_FLAGS) $(INCLUDES) -MD -MP -c $< -o $@

$(BUILD)/rnaseq: $(BUILD)/main.o $(BUILD)/table.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

$(BUILD)/test_table: $(BUILD)/test_table.o $(BUILD)/table.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

verify-vendor:
	sha256sum -c third_party/SHA256SUMS

test: $(BUILD)/test_table
	$(BUILD)/test_table

check: verify-vendor test $(BUILD)/rnaseq
	bash tests/integration/check_native.sh $(BUILD)/rnaseq

sanitize:
	$(MAKE) BUILD=build/sanitize CXXFLAGS='-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined' LDFLAGS='-fsanitize=address,undefined' check

-include $(wildcard $(BUILD)/*.d)
