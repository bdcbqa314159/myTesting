#include "../include/myTesting_bits/textFunctions.hpp"
#include "../include/myTesting_bits/testing.hpp"

using namespace std;

string escapeJavascriptString(const string &in)
{
    int n = in.size();
    stringstream s;
    for (int i = 0; i < n; i++)
    {
        char c = in[i];
        if (c == '\'')
        {
            s << "\\\'";
        }
        else if (c == '\"')
        {
            s << "\\\"";
        }
        else if (c == '\\')
        {
            s << "\\\\";
        }
        else if (c == '\t')
        {
            s << "\\t";
        }
        else if (c == '\n')
        {
            s << "\\n";
        }
        else if (c == '\r')
        {
            s << "\\r";
        }
        else
        {
            s << c;
        }
    }
    return s.str();
}

static void testEscapeJavascriptString()
{
    string in = "\"\'\\\r\n\tNot escaped";
    string out = "\\\"\\\'\\\\\\r\\n\\tNot escaped";
    ASSERT(escapeJavascriptString(in) == out);
}

void testTextFunctions()
{
    TEST(testEscapeJavascriptString);
}