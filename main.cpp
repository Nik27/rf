
#include <cstdint>
#include <iostream>
#include <vector>
#include <cassert>
#include <random>
#include <fstream>
#include <chrono>
#include <thread>
#include <set>
#include <tuple>
#include <array>
#include "util.hpp"

using namespace std;

struct Quasigroup {
    uint64_t n = 0;
    vector <uint64_t> f; // F = {f_1, ..., f_n}
    vector <uint64_t> p; // P = {p_1, ..., p_n}
    vector <uint64_t> s; // Квазигруппа задается латинским квадратом
    vector <uint64_t> r; // Квазигруппа задается G = {g_1, ..., g_n}

    uint64_t p_i(uint64_t x, uint64_t y, uint64_t i) const {
        const uint64_t x_i = ith(x, n, i);
        const uint64_t y_i = ith(y, n, i);
        return ith(p[i], 4, x_i << 1ul | y_i);
    }

    uint64_t f_i(uint64_t x, uint64_t y, uint64_t i) const {
        uint64_t z = 0;
        for (uint64_t j = 0; j < n; j++) {
            z <<= 1ul;
            z |= p_i(x, y, j);
        }
        return ith(f[i], (1ul << n), z);
    }

    uint64_t g_i(uint64_t x, uint64_t y, uint64_t i) const {
        if (r.size() > 0) return ith(r[i], (1ul << 2 * n), x << n | y);
        const uint64_t x_i = ith(x, n, i);
        const uint64_t y_i = ith(y, n, i);
        return x_i ^ y_i ^ f_i(x, y, i);
    }

    uint64_t g(uint64_t x, uint64_t y) const {
        if (s.size() > 0) return s[x * n + y];
        uint64_t z = 0;
        for (uint64_t j = 0; j < n; j++) {
            z = z << 1ul | g_i(x, y, j);
        }
        return z;
    }

    void read_p(istream &is, uint64_t &p) {
        p = 0ul;
        char c;
        for (auto i = 0; i < 4; i++) {
            p <<= 1u;
            is >> c;
            p |= c - '0';
        }
    }

    void read_f(istream &is, uint64_t &f) {
        f = 0ul;
        char c;
        uint64_t s = size();
        for (auto i = 0; i < s; i++) {
            f <<= 1u;
            is >> c;
            f |= c - '0';
        }

    }

    void read_family(istream &is) {
        is >> n;
        s.clear();
        p.resize(n);
        for (auto &i : p) read_p(is, i);
        f.resize(n);
        for (auto &i : f) read_f(is, i);
    }

    void read_matrix(istream &is) {
        is >> n;
        s.resize(n * n);
        for (auto &z : s) {
            is >> z;
        }
    }

    void write_matrix(ostream &os) const {
        os << size() << endl;
        const uint64_t s = size();
        for (uint64_t x = 0; x < s; x++) {
            for (uint64_t y = 0; y < s; y++) {
                os << g(x, y) << " ";
            }
            os << endl;
        }
    }

    void write_p(ostream &os, uint64_t p) const {
        os << bitset<4>(p) << endl;
    }

    void write_f(ostream &os, uint64_t f) const {
        const uint64_t s = size();
        os << bitset<64>(f).to_string().substr(64 - s, 64) << endl;
    }

    void write_ps(ostream &os, bool pretty = false) const {
        for (auto i = 0; i < n; i++) {
            if (pretty) os << "p_" << i + 1 << " = ";
            write_p(os, p[i]);
        }
    }

    void write_fs(ostream &os, bool pretty = false) const {
        for (auto i = 0; i < n; i++) {
            if (pretty) os << "f_" << i + 1 << " = ";
            write_f(os, f[i]);
        }
    }

    void write_family(ostream &os) const {
        os << n << endl;
        write_ps(os, false);
        write_fs(os, false);
    }

    uint64_t size() const {
        return s.empty() ? 1ul << n : n;
    }
};

static bool is_regular(uint64_t family, uint64_t n) {
    assert(n > 0);
    const uint64_t k = 1ul << n;
    if (n < 2) return _is_const(family, k);
    if (n < 3) {
        const uint64_t f_0 = extract(family, n, k, 0);
        const uint64_t f_1 = extract(family, n, k, 1);
        if (_is_const(f_0, k)) return _is_dummy(f_1, n, k, 1);
        if (_is_const(f_1, k)) return _is_dummy(f_0, n, k, 0);
        return false;
    }
    for (uint64_t i = 0; i < n; i++) {
        const uint64_t f = extract(family, n, k, i);
        if (!_is_dummy(f, n, k, i)) return false;
    }

    for (uint64_t a = 0; a < k - 1; a++) {
        for (uint64_t b = a + 1; b < k; b++) {
            bool found = false;
            for (uint64_t i = 0; i < n && !found; i++) {
                if ((a ^ b) & (1ul << (n - 1 - i))) {
                    uint64_t f = extract(family, n, k, i);
                    if (_is_equal(f, k, a, b)) found = true;
                }
            }
            if (!found) return false;
        }
    }
    return true;
}

static void generate_fs(const uint64_t n, vector <uint64_t> &families) {
    assert(1 <= n && n <= 3);
    families.clear();
    const uint64_t k = (1ul << n); // 2^n наборов для n переменных
    const uint64_t m = (~0ul) >> (64u - k * n);

    families.push_back(0);
    for (uint64_t family = 1; family < m; family++) {
        if (is_regular(family, n)) families.push_back(family);
    }
    families.push_back(m);
}

/**
 * Генерируем правильное семейство булевых функций F = {f_1, ..., f_n}
 */
static void generate_f(const uint64_t n, vector <uint64_t> &v) {
    assert(n > 0);
    cout << "Generating " << to_string('F', 'f', n, n) << endl;

    random_device rd;
    mt19937 mt(rd());

    if (n > 3) {
        uniform_int_distribution <uint64_t> d(0, 1);
        v.resize(n);
        for (auto &f : v) f = d(mt) ? ~0ul : 0ul;
    }

    vector <uint64_t> families;
    generate_fs(n, families);

    uniform_int_distribution <uint64_t> distribution(0, families.size() - 1);

    const auto i = distribution(mt);

    cout << i + 1 << " of " << families.size() << endl;

    const uint64_t k = 1ul << n;
    to_vector(families[i], k, n, v);
}

static void generate_p(const uint64_t n, vector <uint64_t> &p) {
    assert(n > 0);
    cout << "Generating " << to_string('P', 'p', n, 2) << endl;

    p.resize(n);
    const uint64_t m = (1ul << 2ul) - 1;

    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution <uint64_t> distribution(0, m);
    for (auto i = 0; i < n; i++) {
        p[i] = distribution(mt);
    }
}

static void print(const Quasigroup &quasigroup, const vector <uint64_t> &subquasigroup, ostream &out = cout) {
    out << subquasigroup.size() << endl;
    for (uint i = 0; i < subquasigroup.size(); i++) {
        for (uint j = 0; j < subquasigroup.size(); j++) {
            out << quasigroup.g(i, j) << " ";
        }
        out << endl;
    }
}

static void f_t(const Quasigroup &quasigroup, uint64_t &sq, uint64_t k, uint64_t l) {
    const uint64_t t = 1ul << quasigroup.g(k, l);
    if (!(sq & t)) sq |= t;
}

static uint64_t find_subquasigroup(const Quasigroup &quasigroup, uint64_t sq, uint64_t su) {
    if (sq == su) return su;
    if (__builtin_popcountll(sq) == quasigroup.size()) return sq;
    uint64_t j = __builtin_ctzll(~(~sq ^ su));
    f_t(quasigroup, sq, j, j);

    uint64_t u = su;
    for (uint64_t t = 0; u; t++) {
        if (u & 1ul) {
            f_t(quasigroup, sq, t, j); // l = j
            f_t(quasigroup, sq, j, t); // k = j
        }
        u >>= 1ul;
    }

    su |= 1ul << j;
    return find_subquasigroup(quasigroup, sq, su);
}

/**
 * Ищет нетривиальную подквазигруппу, содержащую по крайней мере k элементов
 */
static vector <uint64_t> find_subquasigroup(const Quasigroup &quasigroup, uint64_t k) {
    assert(1 <= k && k < quasigroup.size());

    uint64_t sq = (1ul << k) - 1;

    uint64_t subquasigroup = find_subquasigroup(quasigroup, sq, 0);

    auto n = __builtin_popcountll(subquasigroup);
    assert(k <= n && n <= quasigroup.size());
    if (n < quasigroup.size()) {
        vector <uint64_t> s;
        for (uint64_t i = 0; subquasigroup; i++) {
            if (subquasigroup & 1ul) s.push_back(i);
            subquasigroup >>= 1;
        }
        return s;
    }
    return vector<uint64_t>();
}

static uint64_t next(uint64_t subset) {
    const uint64_t t = subset | (subset - 1);
    return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(subset) + 1ul)); // NOLINT
}

static bool find_subquasigroup(const Quasigroup &quasigroup, const uint64_t k, vector <set<uint64_t>> &subquasigroups) {
    const uint64_t m = 1ul << quasigroup.size();
    bool has_no_subquasigroups = true;
    for (uint64_t sq = (1ul << k) - 1; sq < m; sq = next(sq)) {
        uint64_t subquasigroup = find_subquasigroup(quasigroup, sq, 0);
        auto n = __builtin_popcountll(subquasigroup);
        assert(k <= n && n <= quasigroup.size());
        if (n < quasigroup.size()) {
            has_no_subquasigroups = false;
            subquasigroups[n].insert(subquasigroup);
        }
    }
    return has_no_subquasigroups;
}

/**
 * Генерирует квазигруппу размером 2^n x 2^n
 */
static void generate(Quasigroup &quasigroup, const uint64_t n) {
    assert(1 <= n && n <= 6);
    quasigroup.n = n;
    generate_f(n, quasigroup.f);
    quasigroup.write_fs(cout, true);
    cout << endl;
    generate_p(n, quasigroup.p);
    quasigroup.write_ps(cout, true);
    cout << endl;
}

/**
 * Запуск алгоритма для заданной квазигруппы
 */
static void brute_force(const Quasigroup &quasigroup, vector <set<uint64_t>> &subquasigroups,
                        vector <chrono::duration<double>> &time) {
    const uint64_t m = quasigroup.size() / 2;
    for (uint64_t k = 1; k <= m; k++) {
        auto s = chrono::steady_clock::now();
        bool has_no_subquasigroups = find_subquasigroup(quasigroup, k, subquasigroups);
        auto f = chrono::steady_clock::now();
        time[k] += f - s;
        if (has_no_subquasigroups) break;
    }
}

/**
 * Перебор всех квазигрупп порожденных правильными семействами булевых функций и выполнение алгоритма на каждой из них
 */
static void brute_force(uint64_t n, bool verbose) {
    assert(1 <= n && n <= 3);
    if (verbose) cout << "n = " << n << endl;

    vector <uint64_t> families;
    generate_fs(n, families);

    if (verbose) cout << "Count of families: " << families.size() << endl;

    const uint64_t k = 1ul << n;
    const uint64_t m = 1ul << 4 * n;

    set <vector<uint64_t>> quasigroups;

    Quasigroup quasigroup;
    quasigroup.n = n;
    for (auto &family : families) {
        to_vector(family, k, n, quasigroup.f);
        for (uint64_t p = 0; p < m; p++) {
            to_vector(p, 4, n, quasigroup.p);
            vector <uint64_t> r(3);
            for (uint64_t x = 0; x < k; x++) {
                for (uint64_t y = 0; y < k; y++) {
                    for (uint64_t i = 0; i < n; i++) {
                        r[i] <<= 1u;
                        r[i] |= quasigroup.g_i(x, y, i);
                    }
                }
            }
            quasigroups.insert(r);
        }
    }

    cout << quasigroups.size() << endl;

    quasigroup.f.clear();
    quasigroup.p.clear();

    vector <uint64_t> results(k);
    vector <chrono::duration<double>> time(k);
    vector <uint64_t> count_a(k);
    vector <uint64_t> count_b(k);


    auto s = chrono::steady_clock::now();

    for (auto &r : quasigroups) {
        quasigroup.r = r;
        vector <set<uint64_t>> subquasigroups(k);
        brute_force(quasigroup, subquasigroups, time);

        for (int i = 1; i <= k / 2; i++) {
            if (!subquasigroups[i].empty()) {
                count_a[i]++;
                results[i] += subquasigroups[i].size();
            }
        }
        for (int i = 1; i <= k / 2; i++) {
            bool has_no_subquasigroups = true;
            for (int j = i; j <= k / 2; j++) {
                if (!subquasigroups[j].empty()) {
                    has_no_subquasigroups = false;
                    break;
                }
            }
            if (has_no_subquasigroups) count_b[i]++;
        }
    }

    auto f = chrono::steady_clock::now();

    cout << endl << "Elapsed time: " << chrono::duration_cast<chrono::milliseconds>(f - s).count() << "ms" << endl;

    for (auto i = 1; i <= k / 2; i++) {
        cout << "k = " << i << "\t" << results[i] << "\t" << count_a[i] << "\t" << count_b[i] << "\t"
             << chrono::duration_cast<chrono::milliseconds>(time[i]).count() << "ms." << endl;
    }

}

int main(int argc, char *argv[]) {
    /**
     * n - порядок правильного семейства булевых функций, задающего квазигруппу порядка 2^n, используется при генерации
     *     по умолчанию n = 2
     * k - минимальный порядок искомой нетривиальной подквазигруппы
     *     по умолчанию k = 1
     */
    uint64_t n, k;
    string family; // Название входного файла, содержащего квазигруппу заданную правильным семейством булевых функций
    string matrix; // Название входного файла, содержащего квазигруппу заданную матрицей
    string output; // Суффикс выходных файлов, по умолчанию равен 'n'
    bool verbose;  // Вывод входных матриц и отладочной информации на экран
    bool all; // Осуществить полный перебор, учитывается только аргумент -n и -v
    uint64_t m = 0ul;
    uint64_t M = 0ul;
    random_device rd;
    int64_t seed = rd();

    n = 2u;
    k = 1u;

    for (auto i = 1; i < argc; i++) {
        string arg(argv[i]);
        if (arg == "-a") {
            all = true;
        } else if (arg == "-o") {
            output = string(argv[++i]);
        } else if (arg == "-m") {
            matrix = string(argv[++i]);
        } else if (arg == "-f") {
            family = string(argv[++i]);
        } else if (arg == "-n") {
            n = (uint64_t) stoi(argv[++i]);
        } else if (arg == "-k") {
            k = (uint64_t) stoi(argv[++i]);
        } else if (arg == "-v") {
            verbose = true;
        } else if (arg == "-g") {
            m = (uint64_t) stoi(argv[++i]);
            if (!M) M = ~0ul;
        } else if (arg == "-G") {
            M = (uint64_t) stoi(argv[++i]);
            if (!m) m = ~0ul;
        } else if (arg == "-s") {
            seed = stoi(argv[++i]);
        }
    }


    if (m > 0 && M > 0) {
        assert(1 <= n && n <= 4);
        const uint64_t k = 1ul << n;
        const uint64_t b = ~0ul >> (64ul - k);
        cout << "s = " << seed << endl;
        cout << "g = " << (int64_t) m << " G = " << (int64_t) M << endl;
        mt19937_64 mt(seed);
        uniform_int_distribution <uint64_t> distribution(0, b);

        uint64_t i, j;
        for (i = 0, j = 0; i < m && j < M; j++) {
            uint64_t family = 0ul;
            for (auto a = 0; a < n; a++) {
                family <<= k;
                family |= distribution(mt);
            }
            if (is_regular(family, n)) {
                i++;
                for (auto a = 0; a < n; a++) {
                    const uint64_t f = extract(family, n, k, a);
                    cout << bitset<64>(f).to_string().substr(64 - k, 64) << " ";
                }
                cout << endl;
            }
        }
        cout << "g = " << i << " G = " << j << endl;
        exit(0);
    }

    if (all) {
        // Осуществляем полный перебор по всем значениям k, 1 <= k <= n
        // Игнорируем все аргументы, кроме v и n - порядка правильного семейства функций
        brute_force(n, verbose);
        exit(0);
    }

    Quasigroup quasigroup;
    if (family.length() > 0) {
        ifstream family_input(family);
        quasigroup.read_family(family_input);
        family_input.close();

        if (verbose) quasigroup.write_family(cout);
    } else if (matrix.length() > 0) {
        ifstream matrix_input(matrix);
        quasigroup.read_matrix(matrix_input);
        matrix_input.close();

        if (verbose) quasigroup.write_matrix(cout);
    } else {
        generate(quasigroup, n);

        if (output.empty()) {
            output = to_string(n);
        }

        ofstream matrix_output("matrix_" + output + ".txt");
        quasigroup.write_matrix(matrix_output);
        matrix_output.close();

        ofstream family_output("family_" + output + ".txt");
        quasigroup.write_family(family_output);
        family_output.close();
    }

    cout << "|Q| = " << quasigroup.size() << ", k = " << k << endl << endl;

    auto s = chrono::steady_clock::now();

    auto subquasigroup = find_subquasigroup(quasigroup, k);

    auto f = chrono::steady_clock::now();

    if (subquasigroup.empty()) {
        cout << "No subquasigroup found for k >= " << k << endl;
    } else {
        print(quasigroup, subquasigroup);
    }

    cout << endl << "Elapsed time: " << chrono::duration_cast<chrono::milliseconds>(f - s).count() << "ms." << endl;
    return 0;
}
