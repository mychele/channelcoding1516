#include <iostream>
#include <vector>
#include <bitset>

/**
 * A BitIo object stores a bitset into a vector of bytes and writes/reads to file
 * Inspired to http://codereview.stackexchange.com/questions/67012/easy-bitset-i-o/119941#119941,
 * with some corrections for N not a multiple of 8
 */

template <std::size_t N>
class BitIo
{
public:

    void push_back(const std::bitset<N>& bs)
    {
        std::vector<Byte> result(NUM_BYTES_PER_BITSET);
        for (std::size_t bit_index = 0; bit_index < N; bit_index++) {
            result[bit_index / 8] |= (bs[bit_index] << (bit_index % 8));
        }
        m_bytes.insert(std::end(m_bytes), std::begin(result), std::end(result));
        m_num_bytes += NUM_BYTES_PER_BITSET;
    }

    std::bitset<N> pop_front()
    {
        std::bitset<N> result;
        for (std::size_t bit_index = 0; bit_index < N; bit_index++) {
            result[bit_index] = ((m_bytes[(bit_index / 8) + m_offset] >> (bit_index % 8)) % 2);
        }
        m_offset += NUM_BYTES_PER_BITSET;
        m_num_bytes -= NUM_BYTES_PER_BITSET;
        return result;
    }

    bool isEmpty() const
    {
        return m_num_bytes < NUM_BYTES_PER_BITSET;
    }

    void clear() 
    {
        m_bytes.clear();
        m_num_bytes = 0;
    }

    std::size_t size() const
    {
        return m_num_bytes;
    }

private:

    using Byte = unsigned char;
    static constexpr std::size_t NUM_BYTES_PER_BITSET = (N + 7) / 8; // it is as ceil((double)N/8)

    template <std::size_t T>
    friend std::ostream& operator<<(std::ostream& os, const BitIo<T>& bio);
    template <std::size_t T>
    friend std::istream& operator>>(std::istream& is, BitIo<T>& bio);

    std::istream& read_file(std::istream& is)
    {
        m_bytes.clear();

        std::streampos current_pos, file_size;
        current_pos = is.tellg();
        is.seekg(0, std::ios::end);
        file_size = is.tellg() - current_pos;
        is.seekg(current_pos, std::ios::beg);

        m_bytes.resize(file_size);
        is.read(reinterpret_cast<char *>(&m_bytes[0]), file_size);

        m_num_bytes += file_size;

        return is;
    }

    std::vector<Byte> m_bytes;
    std::size_t m_offset = 0;
    std::size_t m_num_bytes = 0;
};

template <std::size_t N>
std::ostream& operator<<(std::ostream& os, const BitIo<N>& bio)
{
    for (const auto& byte : bio.m_bytes) {
        os << byte;
    }
    return os;
}

template <std::size_t N>
std::istream& operator>>(std::istream& is, BitIo<N>& bio)
{
    if(!is) {
        is.setstate(std::ios::failbit);
    }
    bio.read_file(is);
    return is;
}

