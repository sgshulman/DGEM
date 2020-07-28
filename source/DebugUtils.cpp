#include "DebugUtils.hpp"
#include <stdexcept>

void dataAssert(const char* description, const char* file, int line)
{
    throw std::invalid_argument(std::string(file) + ":" + std::to_string(line) + ": " + description);
}
