#!/bin/sh

# Format all adapter C and C++ files with clang-format
find ./ \( -iname "*.h" -o -iname "*.hpp" -o -iname "*.c" -o -iname "*.cpp" \) -exec clang-format -i {} \;
