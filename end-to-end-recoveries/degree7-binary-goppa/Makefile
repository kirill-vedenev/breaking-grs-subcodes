CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -pthread
PYTHON ?= python3

.PHONY: cpu gpu check
cpu: build/extract-jets build/libseries.so build/dual-profile
	$(MAKE) -C holdout cpu

gpu:
	$(MAKE) -C holdout gpu

build:
	mkdir -p build

build/extract-jets: recovery/extract_jets.cpp | build
	$(CXX) $(CXXFLAGS) $< -o $@

build/libseries.so: recovery/series.cpp | build
	$(CXX) $(CXXFLAGS) -shared -fPIC $< -o $@

build/dual-profile: tests/dual_profile.cpp | build
	$(CXX) $(CXXFLAGS) $< -o $@

check:
	$(PYTHON) -B scripts/verify.py
