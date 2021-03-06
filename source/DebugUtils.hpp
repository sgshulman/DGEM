#ifndef DEBUG_UTILS_HPP_
#define DEBUG_UTILS_HPP_

#include <string>

void dataAssert(const char* description, const char* file, int line);
void dataAssert(const std::string&& description, const char* file, int line);

#ifndef DATA_ASSERT
#define DATA_ASSERT(EX, DESCRIPTION) (void)((EX) || (dataAssert(DESCRIPTION, __FILE__, __LINE__), false))
#endif

#endif
