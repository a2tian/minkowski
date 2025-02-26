# Compiler and flags
CXX := g++
CXXFLAGS := -g -O2 -std=c++11 -pthread -march=native -fopenmp -Wall -Wextra -I./cc
LDFLAGS := -lstdc++fs -lntl -lgmp -lm -fopenmp

# Directories
SRCDIR := ./cc
BUILDDIR := ./build_cc
TARGET := test

# Source files and object files
SRC := $(wildcard $(SRCDIR)/*.cc) test.cc
SRC := $(filter-out ./cc/bindings.cc, $(SRC))
OBJ := $(SRC:%.cc=$(BUILDDIR)/%.o)

# Default target
all: $(TARGET)

# Link the target
$(TARGET): $(OBJ)
		$(CXX) $(OBJ) -o $@ $(LDFLAGS)

# Compile source files into object files
$(BUILDDIR)/%.o: %.cc | $(BUILDDIR)
		@mkdir -p $(dir $@)
		$(CXX) $(CXXFLAGS) -c $< -o $@

# Create the build directory if it doesn't exist
$(BUILDDIR):
		@mkdir -p $(BUILDDIR)

# Clean up build artifacts
clean:
		rm -rf $(BUILDDIR) $(TARGET)

.PHONY: all clean