#pragma once
#include "utils.h"
#include "clifford22.h"

struct matrix44
{
    static constexpr size_t ndim = 4;
    using thistype = matrix44;
    using mat_type = std::array<std::complex<double>, ndim* ndim>;
    mat_type mat;
    matrix44()
    {
        mat.fill(0);
    }
    matrix44(const std::array<std::complex<double>, ndim* ndim> &mat_)
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

    template<typename T>
    inline auto operator*=(const T& m) const
    {
        for (int i = 0; i < mat.size(); ++i)
        {
            mat[i] *= m;
        }
        return (*this);
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

    inline matrix44 normalize() const
    {
        if (!is_close_to_zero(mat[0])) {
            matrix44 ret = (*this) * std::conj(mat[0]) * (1.0 / std::abs(mat[0]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
        else if (!is_close_to_zero(mat[1])) {
            matrix44 ret = (*this) * std::conj(mat[1]) * (1.0 / std::abs(mat[1]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
        else if (!is_close_to_zero(mat[2])) {

            matrix44 ret = (*this) * std::conj(mat[2]) * (1.0 / std::abs(mat[2]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
        else if (!is_close_to_zero(mat[3])) {

            matrix44 ret = (*this) * std::conj(mat[3]) * (1.0 / std::abs(mat[3]));
            for (auto& element : ret.mat)
                clarify(element);
            return ret;
        }
        else {
            throw std::runtime_error("what matrix?");
        }
    }
};

using matrix44_t = matrix44::mat_type;

inline matrix44 tensor(const matrix22& m1, const matrix22& m2)
{
    matrix44 mat;
    for (int i = 0; i < m1.ndim; ++i) {
        for (int j = 0; j < m1.ndim; ++j) {
            for (int k = 0; k < m2.ndim; ++k) {
                for (int l = 0; l < m2.ndim; ++l) {
                    mat(i * m2.ndim + k, j * m2.ndim + l) = m1(i, j) * m2(k, l);
                }
            }
        }
    }
    return mat;
}

inline matrix44 XI()
{
    return tensor(X(), I()).normalize();
}
inline matrix44 YI()
{
    return tensor(Y(), I()).normalize();
}
inline matrix44 IX()
{
    return tensor(I(), X()).normalize();
}
inline matrix44 IY()
{
    return tensor(I(), Y()).normalize();
}
inline matrix44 YY()
{
    return tensor(H(), I()).normalize();
}
inline matrix44 IH()
{
    return tensor(I(), H()).normalize();
}
inline matrix44 SI()
{
    return tensor(S(), I()).normalize();
}
inline matrix44 IS()
{
    return tensor(I(), S()).normalize();
}
inline matrix44 I_SX()
{
    return tensor(I(), SX()).normalize();
}
inline matrix44 I_SXdag()
{
    return tensor(I(), SXdag()).normalize();
}
inline matrix44 I_SY()
{
    return tensor(I(), SY()).normalize();
}
inline matrix44 I_SYdag()
{
    return tensor(I(), SYdag()).normalize();
}
inline matrix44 SX_I()
{
    return tensor(SX(), I()).normalize();
}
inline matrix44 SXdag_I()
{
    return tensor(SXdag(), I()).normalize();
}
inline matrix44 SY_I()
{
    return tensor(SY(), I()).normalize();
}
inline matrix44 SYdag_I()
{
    return tensor(SYdag(), I()).normalize();
}

inline matrix44 CZ()
{
    matrix44 cz;
    cz(0, 0) = 1;
    cz(1, 1) = 1;
    cz(2, 2) = 1;
    cz(3, 3) = -1;
    return cz.normalize();
}

inline matrix44 Identity44()
{
    matrix44 id;
    id(0, 0) = 1;
    id(1, 1) = 1;
    id(2, 2) = 1;
    id(3, 3) = 1;
    return id.normalize();
}

enum Generator44Enum : int
{
    Generator_I44,
    Generator_XI,
    Generator_IX,
    Generator_YI,
    Generator_IY,
    Generator_SX_I,
    Generator_SXdag_I,
    Generator_I_SX,
    Generator_I_SXdag,
    Generator_SY_I,
    Generator_SYdag_I,
    Generator_I_SY,
    Generator_I_SYdag,
    Generator_CZ,
};

inline std::map<matrix44, std::vector<int>> initialize_clifford44() {
    std::map<matrix44, std::vector<int>> group;
    auto id = Identity44().normalize();
    auto xi = XI().normalize();
    auto ix = IX().normalize();
    auto yi = YI().normalize();
    auto iy = IY().normalize();
    auto sxi = SX_I().normalize();
    auto sxdagi = SXdag_I().normalize();
    auto isx = I_SX().normalize();
    auto isxdag = I_SXdag().normalize();
    auto syi = SY_I().normalize();
    auto sydagi = SYdag_I().normalize();
    auto isy = I_SY().normalize();
    auto isydag = I_SYdag().normalize();
    auto cz = CZ().normalize();
    group.insert({ id, {Generator_I44} });
    group.insert({ xi, {Generator_XI} });
    group.insert({ ix, {Generator_IX} });
    group.insert({ yi, {Generator_YI} });
    group.insert({ iy, {Generator_IY} });
    group.insert({ sxi, {Generator_SX_I} });
    group.insert({ sxdagi, {Generator_SXdag_I} });
    group.insert({ isx, {Generator_I_SX} });
    group.insert({ isxdag, {Generator_I_SXdag} });
    group.insert({ syi, {Generator_SY_I} });
    group.insert({ sydagi, {Generator_SYdag_I} });
    group.insert({ isy, {Generator_I_SY} });
    group.insert({ isydag, {Generator_I_SYdag} });
    group.insert({ cz, {Generator_CZ} });

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
        std::map<matrix44, std::vector<int>> new_elements;

        for (auto [elem, generator_list] : group) {
            MAKE_NEW_ELEMENT(xi, Generator_XI);

            MAKE_NEW_ELEMENT(yi, Generator_YI);

            MAKE_NEW_ELEMENT(ix, Generator_IX);

            MAKE_NEW_ELEMENT(iy, Generator_IY);

            MAKE_NEW_ELEMENT(cz, Generator_CZ);

            MAKE_NEW_ELEMENT(sxi, Generator_SX_I);

            MAKE_NEW_ELEMENT(sxdagi, Generator_SXdag_I);

            MAKE_NEW_ELEMENT(syi, Generator_SY_I);

            MAKE_NEW_ELEMENT(sydagi, Generator_SYdag_I);

            MAKE_NEW_ELEMENT(isx, Generator_I_SX);

            MAKE_NEW_ELEMENT(isxdag, Generator_I_SXdag);

            MAKE_NEW_ELEMENT(isy, Generator_I_SY);

            MAKE_NEW_ELEMENT(isydag, Generator_I_SYdag);
        }
        group.insert(new_elements.begin(), new_elements.end());
        std::cout << "Group size = " << group.size() << "\n";
    } while (added_new);
    return group;
}

struct KeyValueClifford44
{
    matrix44_t arr;
    std::array<int, 10> buffer;
    int size;
};

inline std::vector<KeyValueClifford44> to_serializable_data(const std::map<matrix44, std::vector<int>>& group) {
    std::vector<KeyValueClifford44> m(group.size());
    int i = 0;
    for (auto&& [key, value] : group)
    {
        KeyValueClifford44& object = m[i];
        object.arr = key.mat;
        object.buffer.fill(-1);
        if (value.size() > 10)
            throw std::runtime_error("more than 10");
        for (int i = 0; i < value.size(); ++i)
        {
            object.buffer[i] = value[i];
        }
        object.size = value.size();
        ++i;
    }
    return m;
}

inline std::vector<int> generate_clifford44_multiplication_table(
    const std::map<matrix44, std::vector<int>>& group,
    const std::vector<KeyValueClifford44>& group_list) {

    const int N = group.size();
    std::vector<int> multiplication_table(N * N, 0);

    // 填充乘法表
    std::cout << "Generating table: " << std::endl;
    int i = 0;
    int gridsize = 288;
    int blocksize = N / 288;
#pragma omp parallel for
    for (int block_id = 0; block_id < gridsize; ++block_id) {
        for (int thread_id = 0; thread_id < blocksize; ++thread_id)
        {
            int i = block_id * blocksize + thread_id;
            auto group1 = matrix44(group_list[i].arr);
            // std::cout << "i = " << i << std::endl;
            for (int j = 0; j < N; ++j) {
                auto& group2 = group_list[j].arr;
                auto product = (group1 * group2).normalize();
                auto iter = group.find(product);
                int result = std::distance(group.begin(), iter);
                if (!(matrix44(group_list[result].arr).normalize() == product))
                {
                    throw std::runtime_error("Bad generation.");
                }
                multiplication_table[i * N + j] = result;
            }
        }
        std::cout << "block_id = " << block_id << " finished.\n";
    }
    //for (const auto& row : multiplication_table) {
    //    for (const auto& element : row) {
    //        std::cout << std::setw(4) << element;  // 设置宽度为4，使得输出整齐
    //    }
    //    std::cout << std::endl;  // 每打印完一行后换行
    //}
    return multiplication_table;
}

inline std::vector<int> generate_clifford44_inverse_table(const std::vector<int> &multiplication_table, 
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

inline std::tuple<std::vector<int>, int> rb44_with_correctness_check(
    const std::vector<int> &table, const std::vector<int> &inverse, int clifford_depth,
    std::vector<KeyValueClifford44> &serializable_group_data
    )
{
    int N = 11520;
    std::default_random_engine eng(10086);
    std::uniform_int_distribution<int> ud(0, N - 1);
    std::vector<int> sequence;
    int current = ud(eng);
    sequence.resize(clifford_depth);
    sequence[0] = current;
    matrix44 t = serializable_group_data[current].arr;
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
    matrix44 inv_mat = serializable_group_data[inv].arr;
    matrix44 expect_to_be_i = (inv_mat * t).normalize();

    if (expect_to_be_i == Identity44())
        return { sequence, inv };
    else
        throw std::runtime_error("bad test.");
}


inline std::tuple<std::vector<int>, int> rb44(
    const std::vector<int>& table, const std::vector<int>& inverse, int clifford_depth,
    std::vector<KeyValueClifford44>& serializable_group_data
)
{
    int N = 11520;
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

inline int generate_readable_group_data44()
{
    // for debug
    auto group = initialize_clifford44();
    auto serializable_group_data = to_serializable_data(group);
    FILE* fp;
    fopen_s(&fp, "readable_group44.txt", "w");

    if (!fp)
    {
        std::cout << "File not found." << std::endl;
        return -1;
    }

    for (int i = 0; i < serializable_group_data.size(); ++i)
    {
        auto& data = serializable_group_data[i];
        std::vector<int> vec(data.size);
        for (int j = 0; j < vec.size(); ++j)
        {
            vec[j] = data.buffer[j];
        }
        fprintf_s(fp, "%5d : %s\n", i, vec2str(vec).c_str());
    }
    fclose(fp);

    return 0;
}

int generate_table44()
{
    auto group = initialize_clifford44();

    auto serializable_group_data = to_serializable_data(group);
    int pos_of_identity = std::distance(group.begin(), group.find(Identity44()));
    auto table = generate_clifford44_multiplication_table(group, serializable_group_data);
    int N = group.size();
    auto inverse_table = generate_clifford44_inverse_table(table, N, pos_of_identity);
    FILE* fp;

    fopen_s(&fp, "rb44.dat", "wb");

    if (!fp)
    {
        std::cout << "File not found." << std::endl;
        std::cout << "Generate Failed." << std::endl;
        return -1;
    }

    fwrite(&pos_of_identity, sizeof(int), 1, fp);
    fwrite(table.data(), sizeof(int), table.size(), fp);
    fwrite(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
    fwrite(serializable_group_data.data(), sizeof(KeyValueClifford44),
        serializable_group_data.size(), fp);

    fclose(fp);

    if (!generate_readable_group_data44())
    {
        std::cout << "Generate Failed." << std::endl;

        return -1;
    }

    std::cout << "Generate End" << std::endl;
    return 0;
}


int load_and_generate_inverse_table44()
{
    int N = 11520;
    FILE* fp;
    fopen_s(&fp, "rb44.dat", "rb");
    if (!fp)
    {
        std::cout << "File not found." << std::endl;
        return -1;
    }
    int pos_of_identity;
    std::vector<int> table(N * N);
    std::vector<int> inverse_table(N);
    std::vector<KeyValueClifford44> serializable_group_data(N);
    fread(&pos_of_identity, sizeof(int), 1, fp);
    fread(table.data(), sizeof(int), table.size(), fp);
    // check_unique(table, N);
    std::vector<int> inv_table = generate_clifford44_inverse_table(table, N, pos_of_identity);
    return 0;
}

int testrb44()
{
    int N = 11520;
    FILE* fp;
    fopen_s(&fp, "rb44.dat", "rb");
    if (!fp)
    {
        std::cout << "File not found." << std::endl;
        return -1;
    }
    int pos_of_identity;
    std::vector<int> table(N * N);
    std::vector<int> inverse_table(N);
    std::vector<KeyValueClifford44> serializable_group_data(N);
    fread(&pos_of_identity, sizeof(int), 1, fp);
    fread(table.data(), sizeof(int), table.size(), fp);
    fread(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
    fread(serializable_group_data.data(), sizeof(KeyValueClifford44),
        serializable_group_data.size(), fp);

    std::default_random_engine eng(10086);
    std::uniform_int_distribution<int> ud(0, N - 1);

    int current_mat = ud(eng);
    matrix44 this_matrix = serializable_group_data[current_mat].arr;
    int i = 0;
    while (i++ >= 0)
    {
        std::cout << i << " test pass" << std::endl;
        int new_mat = ud(eng);
        int next_mat = table[current_mat * N + new_mat];
        auto new_matrix = serializable_group_data[new_mat].arr;
        auto next_matrix = serializable_group_data[next_mat].arr;

        if (!((this_matrix * new_matrix).normalize() == matrix44(next_matrix).normalize()))
        {
            throw std::runtime_error("bad computing");
        }
        this_matrix = next_matrix;
        current_mat = next_mat;
    }
    return 0;
}

int load_and_test_rb44()
{
    int N = 11520;
    FILE* fp;
    fopen_s(&fp, "rb44.dat", "rb");
    if (!fp)
    {
        std::cout << "File not found." << std::endl;
        return -1;
    }
    int pos_of_identity;
    std::vector<int> table(N * N);
    std::vector<int> inverse_table(N);
    std::vector<KeyValueClifford44> serializable_group_data(N);
    fread(&pos_of_identity, sizeof(int), 1, fp);
    fread(table.data(), sizeof(int), table.size(), fp);
    fread(inverse_table.data(), sizeof(int), inverse_table.size(), fp);
    fread(serializable_group_data.data(), sizeof(KeyValueClifford44),
        serializable_group_data.size(), fp);

    std::chrono::time_point tp1 = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000; ++i)
    {
        volatile std::tuple<std::vector<int>, int> ret = rb44(table, inverse_table, 1000, serializable_group_data);
    }
    std::chrono::time_point tp2 = std::chrono::steady_clock::now();

    auto duration = tp2 - tp1;
    std::cout << "RB 2 qubit with 1000 clifford depth\n";
    std::cout << "Generate 1000 random configurations duration = " << std::chrono::duration_cast<std::chrono::microseconds>(duration).count() << "us\n";
    return 0;
}
