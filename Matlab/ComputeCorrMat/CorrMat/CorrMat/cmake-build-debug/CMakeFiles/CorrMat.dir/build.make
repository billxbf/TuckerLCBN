# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/billxu/Desktop/CorrMat/CorrMat

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/CorrMat.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CorrMat.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CorrMat.dir/flags.make

CMakeFiles/CorrMat.dir/main.cpp.o: CMakeFiles/CorrMat.dir/flags.make
CMakeFiles/CorrMat.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CorrMat.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CorrMat.dir/main.cpp.o -c /Users/billxu/Desktop/CorrMat/CorrMat/main.cpp

CMakeFiles/CorrMat.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CorrMat.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/billxu/Desktop/CorrMat/CorrMat/main.cpp > CMakeFiles/CorrMat.dir/main.cpp.i

CMakeFiles/CorrMat.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CorrMat.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/billxu/Desktop/CorrMat/CorrMat/main.cpp -o CMakeFiles/CorrMat.dir/main.cpp.s

# Object files for target CorrMat
CorrMat_OBJECTS = \
"CMakeFiles/CorrMat.dir/main.cpp.o"

# External object files for target CorrMat
CorrMat_EXTERNAL_OBJECTS =

CorrMat: CMakeFiles/CorrMat.dir/main.cpp.o
CorrMat: CMakeFiles/CorrMat.dir/build.make
CorrMat: CMakeFiles/CorrMat.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable CorrMat"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CorrMat.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CorrMat.dir/build: CorrMat

.PHONY : CMakeFiles/CorrMat.dir/build

CMakeFiles/CorrMat.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CorrMat.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CorrMat.dir/clean

CMakeFiles/CorrMat.dir/depend:
	cd /Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/billxu/Desktop/CorrMat/CorrMat /Users/billxu/Desktop/CorrMat/CorrMat /Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug /Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug /Users/billxu/Desktop/CorrMat/CorrMat/cmake-build-debug/CMakeFiles/CorrMat.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CorrMat.dir/depend
