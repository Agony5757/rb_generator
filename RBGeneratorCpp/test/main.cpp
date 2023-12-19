#include "utils.h"
#include "clifford22.h"
#include "clifford44.h"

int main44()
{
    generate_table44();
    // generate_readable_group_data44();
    // load_and_generate_inverse_table44();
    testrb44();
    load_and_test_rb44();
    return 0;
}

int main22()
{
    generate_table22();
    // generate_readable_group_data22();
    // load_and_generate_inverse_table22();
    testrb22();
    load_and_test_rb22();
    return 0;
}

int main()
{
    main22();
    // main44();
}