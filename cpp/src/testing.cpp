#include "../include/myTesting_bits/testing.hpp"

static bool debugEnabled = true;

bool isDebugEnabled()
{
    return debugEnabled;
}

void setDebugEnabled(bool enable)
{
    debugEnabled = enable;
}