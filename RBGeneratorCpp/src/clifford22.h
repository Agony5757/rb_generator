#pragma once
#include "utils.h"

constexpr int clifford22_special_operation_count = 7;

struct matrix22
{
    static constexpr size_t ndim = 2;
    using thistype = matrix22;
    using mat_type = std::array<std::complex<double>, ndim* ndim>;
    std::array<std::complex<double>, ndim * ndim> mat;
    matrix22()
    {
        mat.fill(0);
    }
    matrix22(const std::array<std::complex<double>, ndim* ndim>& mat_)
        :mat(mat_)
    {    }

    inline std::complex<double>& operator()(int x, int y)
    {
        return mat[x + y * ndim];
    }
    inline const std::complex<double>& operator()(int x, int y) const
    {
        return mat[x + y * ndim];
    }

    inline thistype operator*(const thistype& other) const
    {
        thistype ret;
        for (int i = 0; i < ndim; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                for (int k = 0; k < ndim; ++k)
                {
                    ret(i, j) += (*this)(i, k) * other(k, j);
                }
            }
        }
        return ret;
    }

    inline thistype operator*(const mat_type& other) const
    {
        thistype ret;
        for (int i = 0; i < ndim; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                for (int k = 0; k < ndim; ++k)
                {
                    ret(i, j) += (*this)(i, k) * other[k + j * ndim];
                }
            }
        }
        return ret;
    }

    inline bool operator==(const thistype& m) const
    {
        for (size_t i = 0; i < mat.size(); ++i)
        {
            if (!is_close(mat[i], m.mat[i]))
                return false;
        }
        return true;
    }

    inline bool operator<(const thistype& m) const
    {
        return mat < m.mat;
    }

    template<typename T>
    inline auto operator*(const T& m) const
    {
        thistype ret;
        for (int i = 0; i < mat.size(); ++i)
        {
            ret.mat[i] = mat[i] * m;
        }
        return ret;
    }

    inline void print() const
    {
        for (int i = 0; i < ndim; ++i)
        {
            for (int j = 0; j < ndim; ++j)
            {
                std::cout << (*this)(i, j) << "\t";
            }
            std::cout << "\n";
        }
    }

    inline matrix22 normalize() const
    {
        if (is_close_to_zero(mat[0])) {
            matrix22 ret = (*this) * std::conj(mat[1]) * (1.0 / std::abs(mat[1]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
        else {
            matrix22 ret = (*this) * std::conj(mat[0]) * (1.0 / std::abs(mat[0]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
    }
};

using matrix22_t = decltype(matrix22::mat);

inline matrix22 I()
{
    matrix22 id;
    id(0, 0) = 1;
    id(1, 1) = 1;
    return id;
}

inline matrix22 X()
{
    matrix22 x;
    x(0, 1) = 1;
    x(1, 0) = 1;
    return x;
}

inline matrix22 Y()
{
    using namespace std::complex_literals;
    matrix22 x;
    x(0, 1) = -1i;
    x(1, 0) = 1i;
    return x;
}

inline matrix22 SX()
{
    matrix22 sx;
    sx(0, 0) = SQRT2;
    sx(0, 1) = { 0, -SQRT2 };
    sx(1, 0) = { 0, -SQRT2 };
    sx(1, 1) = SQRT2;
    return sx;
}

inline matrix22 SXdag()
{
    matrix22 sx;
    sx(0, 0) = SQRT2;
    sx(0, 1) = { 0, SQRT2 };
    sx(1, 0) = { 0, SQRT2 };
    sx(1, 1) = SQRT2;
    return sx;
}

inline matrix22 SY()
{
    matrix22 sy;
    sy(0, 0) = SQRT2;
    sy(0, 1) = -SQRT2;
    sy(1, 0) = SQRT2;
    sy(1, 1) = SQRT2;
    return sy;
}

inline matrix22 SYdag()
{
    matrix22 sy;
    sy(0, 0) = SQRT2;
    sy(0, 1) = SQRT2;
    sy(1, 0) = -SQRT2;
    sy(1, 1) = SQRT2;
    return sy;
}

inline matrix22 S()
{
    using namespace std::complex_literals;
    matrix22 x;
    x(0, 0) = 1;
    x(1, 1) = 1i;
    return x;
}

inline matrix22 H()
{
    matrix22 x;
    x(0, 0) = SQRT2;
    x(0, 1) = SQRT2;
    x(1, 0) = SQRT2;
    x(1, 1) = -SQRT2;
    return x;
}


enum Generator22Enum : int
{
    Generator_I22,
    Generator_X,
    Generator_Y,
    Generator_SX,
    Generator_SY,
    Generator_SXdag,
    Generator_SYdag,
};

inline std::map<matrix22, std::vector<int>> initialize_clifford22() {
    std::map<matrix22, std::vector<int>> group;
    auto id = I().normalize();
    auto x = X().normalize();
    auto y = Y().normalize();
    auto sx = SX().normalize();
    auto sy = SY().normalize();
    auto sxdag = SXdag().normalize();
    auto sydag = SYdag().normalize();
    group.insert({ id, {Generator_I22} });
    group.insert({ x, {Generator_X} });
    group.insert({ y, {Generator_Y} });
    group.insert({ sx, {Generator_SX} });
    group.insert({ sy, {Generator_SY} });
    group.insert({ sxdag, {Generator_SXdag} });
    group.insert({ sydag, {Generator_SYdag} });

#define MAKE_NEW_ELEMENT(generator_mat, generator) {\
auto new_elem = (elem * generator_mat);\
new_elem = new_elem.normalize();\
auto iter_group_find_##new_elem = group.find(new_elem);\
if (iter_group_find_##new_elem == group.end()) {\
    added_new = true;\
    std::vector<int> new_generator_list = generator_list;\
    new_generator_list.push_back(generator);\
    new_elements.insert({ new_elem, new_generator_list});\
}}

    bool added_new;
    do {
        added_new = false;
        std::map<matrix22, std::vector<int>> new_elements;

        for (auto [elem, generator_list] : group) {

            MAKE_NEW_ELEMENT(x, Generator_X);

            MAKE_NEW_ELEMENT(y, Generator_Y);

            MAKE_NEW_ELEMENT(sx, Generator_SX);

            MAKE_NEW_ELEMENT(sy, Generator_SY);

            MAKE_NEW_ELEMENT(sxdag, Generator_SXdag);

            MAKE_NEW_ELEMENT(sydag, Generator_SYdag);
        }
        group.insert(new_elements.begin(), new_elements.end());
        std::cout << "Group size = " << group.size() << "\n";
    } while (added_new);
    return group;
}

struct KeyValueClifford22
{
    matrix22_t arr;
    std::array<int, 3> buffer;
    int size;
};

inline void check_unique(const std::vector<int>& multiplication_table, int N)
{
    for (int i = 0; i < N; ++i)
    {
        std::cout << "i = " << i;
        std::set<int> s;
        for (int j = 0; j < N; ++j)
        {
            int num = multiplication_table[i * N + j];
            if (num < 0 || num >= N)
                throw std::runtime_error("bad number");

            if (s.find(num) != s.end())
                throw std::runtime_error("repeated number");

            s.insert(num);
        }
        if (s.size() != N)
            throw std::runtime_error("repeated number");

        auto iter = s.find(N - 1);
        std::cout << "   Distance = " << std::distance(s.begin(), iter) << std::endl;
    }
}

inline std::vector<KeyValueClifford22> to_serializable_data(const std::map<matrix22, std::vector<int>>& group) {
    std::vector<KeyValueClifford22> m(group.size());
    int i = 0;
    for (auto&& [key, value] : group)
    {
        KeyValueClifford22& object = m[i];
        object.arr = key.mat;
        object.buffer.fill(-1);
        if (value.size() > 3)
            throw std::runtime_error("more than 3");
        for (int i = 0; i < value.size(); ++i)
        {
            object.buffer[i] = value[i];
        }
        object.size = value.size();
        ++i;
    }
    return m;
}

inline std::vector<int> generate_clifford22_multiplication_table(
    const std::map<matrix22, std::vector<int>>& group,
    const std::vector<KeyValueClifford22>& group_list) {

    const int N = group.size();
    std::vector<int> multiplication_table(N * N, 0);

    // ���˷���
    std::cout << "Generating table: " << std::endl;
    int i = 0;
    for (int i = 0; i < N; ++i) {
        auto group1 = matrix22(group_list[i].arr);
        // std::cout << "i = " << i << std::endl;
        for (int j = 0; j < N; ++j) {
            auto& group2 = group_list[j].arr;
            auto product = (group1 * group2).normalize();
            auto iter = group.find(product);
            int result = std::distance(group.begin(), iter);
            if (!(matrix22(group_list[result].arr).normalize() == product))
            {
                throw std::runtime_error("Bad generation.");
            }
            multiplication_table[i * N + j] = result;
        }
        std::cout << "i = " << i << " finished.\n";
    }
    return multiplication_table;
}

inline std::vector<int> generate_clifford22_inverse_table(const std::vector<int>& multiplication_table,
    int N, int pos_of_identity)
{
    std::vector<int> inverse_table(N, 0);
    for (int i = 0; i < N; ++i)
    {
        auto line_begin = multiplication_table.begin() + N * i;
        auto line_end = multiplication_table.begin() + N * (i + 1);
        auto iter = std::find(line_begin, line_end, pos_of_identity);
        if (iter == line_end)
        {
            throw std::runtime_error("Cannot find in line");
        }
        inverse_table[i] = std::distance(line_begin, iter);
    }
    return inverse_table;
}

inline std::tuple<std::vector<int>, int> rb22_with_correctness_check(
    const std::vector<int>& table, const std::vector<int>& inverse, int clifford_depth,
    int N, std::vector<KeyValueClifford22>& serializable_group_data
)
{
    std::default_random_engine eng(10086);
    std::uniform_int_distribution<int> ud(0, N - 1);
    std::vector<int> sequence;
    int current = ud(eng);
    sequence.resize(clifford_depth);
    sequence[0] = current;
    matrix22 t = serializable_group_data[current].arr;
    for (int i = 1; i < clifford_depth; ++i)
    {
        // generate the next
        sequence[i] = ud(eng);
        current = table[current * N + sequence[i]];
        t = t * (serializable_group_data[sequence[i]].arr);
        t = t.normalize();
    }

    // find the inverse
    int inv = inverse[current];
    matrix22 inv_mat = serializable_group_data[inv].arr;
    matrix22 expect_to_be_i = (inv_mat * t).normalize();

    if (expect_to_be_i == I())
        return { sequence, inv };
    else
        throw std::runtime_error("bad test.");
}


inline std::tuple<std::vector<int>, int> rb22(
    const std::vector<int>& table, const std::vector<int>& inverse, int clifford_depth,
    int N, std::vector<KeyValueClifford22>& serializable_group_data
)
{
    std::default_random_engine eng(10086);
    std::uniform_int_distribution<int> ud(0, N - 1);
    std::vector<int> sequence;
    int current = ud(eng);
    sequence.resize(clifford_depth);
    sequence[0] = current;
    for (int i = 1; i < clifford_depth; ++i)
    {
        // generate the next
        sequence[i] = ud(eng);
        current = table[current * N + sequence[i]];
    }

    // find the inverse
    int inv = inverse[current];
    return { sequence, inv };
}

/**************************************************************/
/**************************************************************/
/*********************                *************************/
/********************* Main Functions *************************/
/*********************                *************************/
/**************************************************************/
/**************************************************************/

// inline int generate_readable_group_data22()
// {
//     // for debug
//     auto group = initialize_clifford22();
//     auto serializable_group_data = to_serializable_data(group);
//     FILE* fp;
//     fopen_s(&fp, "readable_group22.txt", "w");

//     if (!fp)
//     {
//         std::cout << "File not found." << std::endl;
//         return -1;
//     }

//     for (int i = 0; i < serializable_group_data.size(); ++i)
//     {
//         auto& data = serializable_group_data[i];
//         std::vector<int> vec(data.size);
//         for (int j = 0; j < vec.size(); ++j)
//         {
//             vec[j] = data.buffer[j];
//         }
//         fprintf_s(fp, "%5d : %s\n", i, vec2str(vec).c_str());
//     }
//     fclose(fp);

//     return 0;
// }

// int generate_table22()
// {
//     auto group = initialize_clifford22();
//     auto id = I().normalize();
//     auto x = X().normalize();
//     auto y = Y().normalize();
//     auto sx = SX().normalize();
//     auto sy = SY().normalize();
//     auto sxdag = SXdag().normalize();
//     auto sydag = SYdag().normalize();
//     auto serializable_group_data = to_serializable_data(group);
//     std::vector<int> special_operation_table;
//     special_operation_table.push_back(std::distance(group.begin(), group.find(id)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(x)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(y)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(sx)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(sy)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(sxdag)));
//     special_operation_table.push_back(std::distance(group.begin(), group.find(sydag)));
//     auto table = generate_clifford22_multiplication_table(group, serializable_group_data);
//     int N = group.size();
//     auto inverse_table = generate_clifford22_inverse_table(table, N, special_operation_table[0]);
//     FILE* fp;

//     fopen_s(&fp, "rb22.dat", "wb");

//     if (!fp)
//     {
//         std::cout << "File not found." << std::endl;
//         std::cout << "Generate Failed." << std::endl;
//         return -1;
//     }

//     fwrite(special_operation_table.data(), sizeof(int), special_operation_table.size(), fp);
//     fwrite(table.data(), sizeof(int), table.size(), fp);
//     fwrite(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
//     fwrite(serializable_group_data.data(), sizeof(KeyValueClifford22),
//         serializable_group_data.size(), fp);

//     fclose(fp);

//     if (generate_readable_group_data22())
//     {
//         std::cout << "Generate Failed." << std::endl;

//         return -1;
//     }

//     std::cout << "Generate End" << std::endl;
//     return 0;
// }


// int load_and_generate_inverse_table22()
// {
//     int N = 24;
//     FILE* fp;
//     fopen_s(&fp, "rb22.dat", "rb");
//     if (!fp)
//     {
//         std::cout << "File not found." << std::endl;
//         return -1;
//     }
//     std::vector<int> special_operations(clifford22_special_operation_count);
//     std::vector<int> table(N * N);
//     std::vector<int> inverse_table(N);
//     std::vector<KeyValueClifford22> serializable_group_data(N);
//     fread(special_operations.data(), sizeof(int), special_operations.size(), fp);
//     fread(table.data(), sizeof(int), table.size(), fp);
//     // check_unique(table, N);
//     std::vector<int> inv_table = generate_clifford22_inverse_table(table, N, special_operations[0]);
//     return 0;
// }

// int testrb22()
// {
//     int N = 24;
//     FILE* fp;
//     fopen_s(&fp, "rb22.dat", "rb");
//     if (!fp)
//     {
//         std::cout << "File not found." << std::endl;
//         return -1;
//     }
//     std::vector<int> special_operations(clifford22_special_operation_count);
//     std::vector<int> table(N * N);
//     std::vector<int> inverse_table(N);
//     std::vector<KeyValueClifford22> serializable_group_data(N);
//     fread(special_operations.data(), sizeof(int), special_operations.size(), fp);
//     fread(table.data(), sizeof(int), table.size(), fp);
//     fread(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
//     fread(serializable_group_data.data(), sizeof(KeyValueClifford22),
//         serializable_group_data.size(), fp);

//     std::default_random_engine eng(10086);
//     std::uniform_int_distribution<int> ud(0, N - 1);

//     int current_mat = ud(eng);
//     matrix22 this_matrix = serializable_group_data[current_mat].arr;
//     int i = 1000;
//     while (i --> 0)
//     {
//         int new_mat = ud(eng);
//         int next_mat = table[current_mat * N + new_mat];
//         auto new_matrix = serializable_group_data[new_mat].arr;
//         auto next_matrix = serializable_group_data[next_mat].arr;

//         if (!((this_matrix * new_matrix).normalize() == matrix22(next_matrix).normalize()))
//         {
//             throw std::runtime_error("bad computing");
//         }
//         this_matrix = next_matrix;
//         current_mat = next_mat;
//     }
//     std::cout << "test passed" << std::endl;
//     return 0;
// }

// inline bool rb22_checker(const std::vector<int> &sequence, const std::vector<int> &inv_sequence)
// {
//     matrix22 m = I();

//     for (auto gate : sequence)
//     {
//         switch (gate)
//         {    
//         case Generator_I22:
//             m = m * I(); break;
//         case Generator_X:
//             m = m * X(); break;
//         case Generator_Y:
//             m = m * Y(); break;
//         case Generator_SX:
//             m = m * SX(); break;
//         case Generator_SY:
//             m = m * SY(); break;
//         case Generator_SXdag:
//             m = m * SXdag(); break;
//         case Generator_SYdag:
//             m = m * SYdag(); break;
//         default:
//             throw std::runtime_error("Bad gate in sequence");
//         }
//     }
//     for (auto gate : inv_sequence)
//     {
//         switch (gate)
//         {
//         case Generator_I22:
//             m = m * I(); break;
//         case Generator_X:
//             m = m * X(); break;
//         case Generator_Y:
//             m = m * Y(); break;
//         case Generator_SX:
//             m = m * SX(); break;
//         case Generator_SY:
//             m = m * SY(); break;
//         case Generator_SXdag:
//             m = m * SXdag(); break;
//         case Generator_SYdag:
//             m = m * SYdag(); break;
//         default:
//             throw std::runtime_error("Bad gate in inv_sequence");
//         }
//     }
//     if (m.normalize() == I())
//         return true;
//     else
//         return false;
// }


struct RB22
{
    int n_special_operations = 7;
    int N = 24;
    std::vector<int> special_operations;
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
        FILE* fp = fopen(filename.c_str(), "rb");
        if (!fp)
        {
            throw std::runtime_error("File not found.");
        }
        table = std::vector<int>(N * N);
        inverse_table = std::vector<int>(N);
        serializable_group_data = std::vector<KeyValueClifford22>(N);
        special_operations = std::vector<int>(clifford22_special_operation_count);
        fread(special_operations.data(), sizeof(int), n_special_operations, fp);
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
        return get_group_elements()[x].arr;
    }

    inline const auto& get_generator(int x) const
    {
        return get_group_elements()[x].buffer;
    }

    inline int get_generator_size(int x) const
    {
        return get_group_elements()[x].size;
    }

    inline std::tuple<std::vector<int>, std::vector<int>>
        get_full_sequence_and_inverse_sequence(const std::vector<int>& input_sequence) const
    {
        check_loaded();
        std::vector<int> full_sequence;
        int multiplies = -1;
        for (int input : input_sequence)
        {
            if (multiplies == -1)
                multiplies = input;
            else
                multiplies = table[multiplies * N + input];

            const auto& generator = serializable_group_data[input].buffer;
            int generator_size = serializable_group_data[input].size;
            full_sequence.insert(full_sequence.end(), generator.begin(), generator.begin() + generator_size);
        }
        int inverse = inverse_table[multiplies];
        const auto& inverse_generator = serializable_group_data[inverse].buffer;
        int inverse_generator_size = serializable_group_data[inverse].size;
        std::vector<int> inverse_sequence;
        inverse_sequence.insert(inverse_sequence.end(), inverse_generator.begin(), inverse_generator.begin() + inverse_generator_size);
        return { full_sequence, inverse_sequence };
    }

    inline std::vector<std::string> get_special_operations_str() const
    {
        return { "I", "X", "Y", "X/2", "Y/2", "-X/2", "-Y/2" };
    }

    inline const std::vector<int>& get_special_operations() const
    {
        check_loaded();
        return special_operations;
    }

};

int load_and_test_rb22(int n_episode = 1000, int clifford_length = 1000)
{
    RB22 rb22;
    rb22.load_from_file("rb22.dat");

    std::chrono::time_point tp1 = std::chrono::steady_clock::now();
    std::uniform_int_distribution<int> ud(0, rb22.N - 1);
    std::default_random_engine eng(10085);
    std::vector<int> random_sequence(clifford_length);
    for (int i = 0; i < n_episode; ++i)
    {
        for (int j = 0; j < clifford_length; ++j)
            random_sequence[j] = ud(eng);
        auto&& [seq, inv] = rb22.get_full_sequence_and_inverse_sequence(random_sequence);
        if (!rb22_checker(seq, inv))
            throw std::runtime_error("Checker not passed.");
    }
    std::chrono::time_point tp2 = std::chrono::steady_clock::now();

    auto duration = tp2 - tp1;
    std::cout << "RB 1 qubit with 1000 clifford depth\n";
    std::cout << "Generate 1000 random configurations duration = " << std::chrono::duration_cast<std::chrono::microseconds>(duration).count() << "us\n";
    return 0;
}