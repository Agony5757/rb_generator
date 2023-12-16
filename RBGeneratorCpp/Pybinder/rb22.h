
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"
#include "pybind11/operators.h"
#include "pybind11/numpy.h"
#include "clifford22.h"

namespace py = pybind11;
using namespace pybind11::literals;

struct RB22
{
    int N = 24;
    int pos_of_identity;
    std::vector<int> table;
    std::vector<int> inverse_table;
    std::vector<KeyValueClifford22> serializable_group_data;
    bool loaded = false;

    RB22() {}

    inline void check_loaded() const
    {
        if (!loaded) 
            throw std::runtime_error("Load data before using RB classes!");
    }

    inline bool load_from_file(const std::string& filename)
    {
        FILE* fp = fopen("rb22.dat", "rb");
        if (!fp)
        {
            throw std::runtime_error("File not found.");
        }
        table = std::vector<int>(N * N);
        inverse_table = std::vector<int>(N);
        serializable_group_data = std::vector<KeyValueClifford22>(N);
        fread(&pos_of_identity, sizeof(int), 1, fp);
        fread(table.data(), sizeof(int), table.size(), fp);
        fread(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
        fread(serializable_group_data.data(), sizeof(KeyValueClifford22),
            serializable_group_data.size(), fp);
        fclose(fp);
        loaded = true;
        return true;
    }

    int get_multiplication_table_elem(int x, int y) const
    {
        check_loaded();
        return table[x * N + y];
    }

    int get_inverse(int x) const
    {
        check_loaded();
        return inverse_table[x];
    }

    inline const std::vector<KeyValueClifford22>& get_group_elements() const
    {
        check_loaded();
        return serializable_group_data;
    }

    inline const std::vector<int>& get_table() const
    {
        check_loaded();
        return table;
    }

    inline const std::vector<int>& get_inverse_table() const
    {
        check_loaded();
        return inverse_table;
    }

    inline const matrix22_t& get_matrix(int x) const
    {
        check_loaded();
        return get_group_elements()[x].arr;
    }

    inline const auto& get_generator(int x) const
    {
        check_loaded();
        return get_group_elements()[x].buffer;
    }

    inline int get_generator_size(int x) const
    {
        check_loaded();
        return get_group_elements()[x].size;
    }

};
