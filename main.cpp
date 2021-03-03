#include <vector>
#include <cstdlib>
#include <cstring>
#include <string>

class BigInt {
    typedef uint32_t T;
    typedef uint64_t GlobT;

    GlobT base;
    bool negative;
    std::vector<T> container;

    void forwardAdd(unsigned curPos, T num) {
        if (num == 0)
            return;
        if (curPos >= container.size())
            container.push_back(0);

        if (T(container[curPos] + num) > T(container[curPos])) {
            container[curPos] += num;
            return;
        }
        container[curPos] += num;
        forwardAdd(curPos + 1, 1);
    }

    void forwardSub(unsigned curPos, T num) {
        if (num == 0)
            return;
        if (curPos >= container.size())
            return;
        if (T(container[curPos]) >= T(num)) {
            container[curPos] -= num;
            return;
        }
        container[curPos] = T(GlobT (base) - GlobT (num) + GlobT(container[curPos]));
        forwardSub(curPos + 1, 1);
    }

    void trailingZerosRm() {
        for (int i = int(size()) - 1; i > 0; --i) {
            if (container[i] != 0)
                return;
            container.pop_back();
        }
        if (container.empty()) {
            container.push_back(0);
            negative = false;
        }
    }

public:
    explicit BigInt(unsigned cSize) :
            container(cSize, 0),
            negative(false),
            base(GlobT(T(-1)) + 1) {}

    BigInt() :
            container(1, 0),
            negative(false),
            base(GlobT(T(-1)) + 1) {}

    BigInt(BigInt const &other) = default;

    BigInt(BigInt &&that) noexcept:
            container(that.container),
            negative(that.negative),
            base(that.base) {
        that.container.clear();
    }

    BigInt &operator=(const BigInt &other) {
        if (this == &other) return *this;
        container = other.container;
        negative = other.negative;
        base = other.base;
        return *this;
    }

    ~BigInt() = default;

    [[nodiscard]] auto begin() {
        return container.begin();
    }

    [[nodiscard]] auto end() {
        return container.end();
    }

    [[nodiscard]] size_t size() const {
        return container.size();
    }

    void zero() {
        negative = false;
        container.clear();
        container.push_back(0);
    }

    BigInt &operator=(int other) {
        zero();
        negative = other < 0;
        container[0] = other;
        return *this;
    }

    BigInt &operator=(long long other) {
        zero();
        negative = other < 0;
        container[0] = other;
        return *this;
    }

    BigInt &operator=(unsigned other) {
        zero();
        negative = false;
        container[0] = other;
        return *this;
    }

    BigInt &operator=(unsigned long long other) {
        zero();
        negative = false;
        container[0] = other;
        return *this;
    }

    int compare(const BigInt &other) const { //0 this == a || -1 this < a || 1 this > a
        if (!negative and other.negative)
            return 1;
        if (!other.negative and negative)
            return -1;
        int check = 1;
        if (negative)
            check = -1;

        if (size() < other.size()) return -1 * check;
        if (size() > other.size()) return check;

        for (size_t i = 0; i < size(); ++i) {
            if (container[i] < other.container[i]) return -1 * check;
            if (container[i] > other.container[i]) return check;
        }
        return 0;
    }

    bool operator<(const BigInt &other) const {
        return compare(other) == -1;
    }

    bool operator>(const BigInt &other) const {
        return compare(other) == 1;
    }

    bool operator==(const BigInt &other) const {
        return compare(other) == 0;
    }

    bool operator!=(const BigInt &other) const {
        return compare(other) != 0;
    }

    bool operator<=(const BigInt &other) const {
        auto res = compare(other);
        return res == 0 || res == -1;
    }

    bool operator>=(const BigInt &other) const {
        auto res = compare(other);
        return res == 0 || res == 1;
    }

    BigInt operator+(const BigInt &other) const {
        auto res = BigInt(*this);
        if (other.negative == res.negative) {
            for (unsigned i = 0; i < other.size(); i++)
                res.forwardAdd(i, other.container[i]);
            res.trailingZerosRm();
            return res;
        }
        return *this - (-other);
    }

    BigInt &operator+=(const BigInt &other) {
        *this = *this + other;
        return *this;
    }

    BigInt operator-(const BigInt &other) const {
        if (negative == other.negative) {
            if (abs() >= other.abs()) {
                BigInt res = *this;
                for (unsigned i = 0; i < other.size(); i++)
                    res.forwardSub(i, other.container[i]);
                res.trailingZerosRm();
                return res;
            }
            return -(other - *this);
        }
        return *this + (-other);
    }

    BigInt &operator-=(const BigInt &other) {
        *this = *this - other;
        return *this;
    }

    BigInt &operator*=(T other) {
        negative = negative xor (other < 0);
        GlobT cary = 0;

        for (int i = 0; i < size(); ++i) {
            GlobT res = cary + GlobT(container[i]) * GlobT(other);
            container[i] = res % base;
            cary = res / base;
        }

        if (cary != 0)
            container.push_back(cary);

        return *this;
    }

    BigInt operator*(T other) {
        BigInt res = *this;
        *this *= other;
        return res;
    }

    BigInt operator-() const {
        BigInt res = *this;
        res.negative = !res.negative;
        return res;
    }

    BigInt &operator+=(unsigned other) {
        forwardAdd(0, other);
        return *this;
    }

    BigInt abs() const {
        BigInt res = *this;
        res.negative = false;
        return res;
    }

    BigInt operator+(unsigned other) const {
        auto retValue = BigInt(*this);
        retValue += other;
        return retValue;
    }

    BigInt& operator-=(int other) {
        auto t = BigInt(1);
        t = other;
        operator-=(t);
        return *this;
    }

    BigInt operator-(int other) {
        auto t = BigInt(*this);
        t -= other;
        return t;
    }

    BigInt &operator=(const std::string &other) {
        zero();
        long ch = 0;
        if (other[0] == '-') {
            ch = 1;
            negative = true;
        }

        for (; ch < other.length(); ch++) {
            operator*=(10);
            operator+=(other[ch] - '0');
        }

        return *this;
    }

    BigInt &shR1() {
        GlobT carry = 0;
        const GlobT hi_bit_set = T(1) << (sizeof(T) * 8 - 1);

        for (auto i = size() - 1; i + 1; --i) {
            const GlobT next_carry = (container[i] & 1) ? hi_bit_set : 0;
            container[i] >>= 1;
            container[i] |= carry;
            carry = next_carry;
        }

        trailingZerosRm();

        return *this;
    }

    void setBitAt(size_t index, bool set = true) {
        size_t widx = index / (sizeof(T) * 8);
        size_t bidx = index % (sizeof(T) * 8);
        if (size() < widx + 1) {
            container.resize(widx + 1);
        }
        if (set) {
            container[widx] |= T(1) << bidx;
        } else {
            container[widx] &= ~(T(1) << bidx);
        }
    }

    BigInt &shL1() {
        GlobT carry = 0;
        const GlobT hi_bit_set = T(1) << (sizeof(T) * 8 - 1);

        for (size_t i = 0; i < size(); ++i) {
            const GlobT next_carry = !!(container[i] & hi_bit_set);
            container[i] <<= 1;
            container[i] |= carry;
            carry = next_carry;
        }

        if (carry) { container.push_back(1); }
        return *this;
    }

    bool isZero() const {
        for (auto c: container){
            if (c!=0)
                return false;
        }
        return true;
    }

    std::pair<BigInt, T> divide(T d) const {
        if(d == 0) {
            // TODO: handle divide by zero
            return {};
        }

        auto result = BigInt(size() + 1);
        result.negative = negative xor (d < 0);
        GlobT cary = 0;

        const GlobT divisor = d;
        for (int i = size() - 1; i >= 0; i--) {
            cary = cary * base + this->container[i];
            result.container[i] = cary / divisor;
            cary %= divisor;
        }

        result.trailingZerosRm();
        return {result, T(cary)};
    }

    BigInt& operator/=(T other) {
        *this = divide(other).first;
        return *this;
    }

    BigInt operator/(T other) const {
        return  divide(other).first;
    }

    T operator%(T other) const{
        return divide(other).second;
    }

    T operator%=(T other){
        *this = divide(other).second;
        return this->container[0];
    }

    std::string to_string() const {
        std::string dec_string;
        auto bi = BigInt(*this);
        T d = 10;

        do {
            const auto next_bi = bi.divide(d);
            const char digit_value = static_cast<char>(next_bi.second);
            dec_string.push_back('0' + digit_value);
            bi = next_bi.first;
        } while(!bi.isZero());
        std::reverse(dec_string.begin(), dec_string.end());
        return dec_string;
    }

    int toInt() const{
        return int(container[0]) * (negative ? -1 : 1);
    }

    unsigned toUnsigned() const{
        return unsigned(container[0]);
    }

    unsigned long long toUnsignedLL() const{
        return (unsigned long long) (container[0]);
    }

    long long toLL() const{
        return (long long) (container[0]) * (negative ? -1 : 1);
    }
};

int main() {
    auto n = BigInt();
    n = 1;
    for (unsigned i = 2 ; i <= 9000; ++i)
        n *= i;

    printf("%s\n", n.to_string().c_str());
    return 0;
}
