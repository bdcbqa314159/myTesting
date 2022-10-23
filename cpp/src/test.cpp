#include "../include/myTesting_bits/test.hpp"
#include "../include/myTesting_bits/standard.hpp"
#include "../include/myTesting_bits/testing.hpp"

using namespace std;

static void message()
{
    cout << endl;
    cout << "myTesting project - my unit test framework." << endl;
}

void testStart()
{
    TEST(message);
}
