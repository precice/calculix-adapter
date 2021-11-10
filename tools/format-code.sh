#!/bin/sh

# Format all adapter C and C++ files with clang-format
find ./adapter/ \( -iname "*.h*" -o -iname "*.c*" \) -exec clang-format -i {} \;
