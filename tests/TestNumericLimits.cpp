#include <iostream>
#include <limits>

int main()
{
    std::cout << std::numeric_limits<double>::min() << std::endl;
    std::cout << std::numeric_limits<double>::lowest() << std::endl;
}