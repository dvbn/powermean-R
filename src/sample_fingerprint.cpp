// SHA-256 and the frozen powermean.sample/1 byte format.
#define R_NO_REMAP
#include <R.h>
#include <R_ext/Memory.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string_view>

#if defined(__APPLE__) && __has_include(<CommonCrypto/CommonDigest.h>)
#include <CommonCrypto/CommonDigest.h>
#define POWERMEAN_COMMONCRYPTO_SHA256
#endif

#include "native_calls.h"

namespace {

using Byte = std::uint8_t;
using Size = std::uint64_t;
using HexDigest = std::array<char, 65>;

static_assert(sizeof(double) == sizeof(Size) && std::numeric_limits<double>::is_iec559,
              "Sample fingerprints require IEEE binary64.");

void add_length(Size& length, std::size_t size) {
    if (size > std::numeric_limits<Size>::max() / 8 - length) {
        Rf_error("The sample fingerprint is outside the SHA-256 length range.");
    }
    length += size;
}

#ifdef POWERMEAN_COMMONCRYPTO_SHA256
class Sha256 {
public:
    Sha256() { CC_SHA256_Init(&context_); }

    void update(const void* source, std::size_t size) {
        add_length(length_, size);
        const auto* bytes = static_cast<const Byte*>(source);
        constexpr std::size_t max_chunk = std::numeric_limits<CC_LONG>::max();
        while (size != 0) {
            const auto take = std::min(size, max_chunk);
            CC_SHA256_Update(&context_, bytes, static_cast<CC_LONG>(take));
            bytes += take;
            size -= take;
        }
    }

    HexDigest finish() {
        std::array<Byte, CC_SHA256_DIGEST_LENGTH> bytes;
        CC_SHA256_Final(bytes.data(), &context_);
        constexpr char digits[] = "0123456789abcdef";
        HexDigest output{};
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            output[2 * i] = digits[bytes[i] >> 4];
            output[2 * i + 1] = digits[bytes[i] & 0xf];
        }
        return output;
    }

private:
    CC_SHA256_CTX context_;
    Size length_ = 0;
};
#else
using Word = std::uint32_t;

// SHA-256 round constants and initial state from FIPS 180-4.
constexpr std::array<Word, 64> rounds{
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

Word rotate_right(Word value, unsigned shift) {
    return (value >> shift) | (value << (32 - shift));
}

Word big_endian_word(const Byte* bytes) {
    return (Word{bytes[0]} << 24) | (Word{bytes[1]} << 16) |
           (Word{bytes[2]} << 8) | Word{bytes[3]};
}

class Sha256 {
public:
    void update(const void* source, std::size_t size) {
        add_length(length_, size);
        if (size == 0) return;
        const auto* bytes = static_cast<const Byte*>(source);
        if (used_ != 0) {
            const auto take = std::min(size, buffer_.size() - used_);
            std::memcpy(buffer_.data() + used_, bytes, take);
            used_ += take;
            bytes += take;
            size -= take;
            if (used_ < buffer_.size()) return;
            compress(buffer_.data());
            used_ = 0;
        }
        while (size >= buffer_.size()) {
            compress(bytes);
            bytes += buffer_.size();
            size -= buffer_.size();
        }
        if (size != 0) std::memcpy(buffer_.data(), bytes, size);
        used_ = size;
    }

    HexDigest finish() {
        buffer_[used_++] = 0x80;
        if (used_ > 56) {
            std::fill(buffer_.begin() + used_, buffer_.end(), 0);
            compress(buffer_.data());
            used_ = 0;
        }
        std::fill(buffer_.begin() + used_, buffer_.begin() + 56, 0);
        const Size bits = length_ * 8;
        for (unsigned i = 0; i < 8; ++i) buffer_[63 - i] = static_cast<Byte>(bits >> (8 * i));
        compress(buffer_.data());

        constexpr char digits[] = "0123456789abcdef";
        HexDigest output{};
        for (std::size_t i = 0; i < state_.size(); ++i) {
            for (unsigned j = 0; j < 8; ++j) {
                output[8 * i + j] = digits[(state_[i] >> (28 - 4 * j)) & 0xf];
            }
        }
        return output;
    }

private:
    void compress(const Byte* bytes) {
        std::array<Word, 64> schedule;
        for (std::size_t i = 0; i < 16; ++i) schedule[i] = big_endian_word(bytes + 4 * i);
        for (std::size_t i = 16; i < schedule.size(); ++i) {
            const Word x = schedule[i - 15];
            const Word y = schedule[i - 2];
            const Word first = rotate_right(x, 7) ^ rotate_right(x, 18) ^ (x >> 3);
            const Word second = rotate_right(y, 17) ^ rotate_right(y, 19) ^ (y >> 10);
            schedule[i] = schedule[i - 16] + first + schedule[i - 7] + second;
        }
        auto [a, b, c, d, e, f, g, h] = state_;
        for (std::size_t i = 0; i < rounds.size(); ++i) {
            const Word choose = (e & f) ^ (~e & g);
            const Word majority = (a & b) ^ (a & c) ^ (b & c);
            const Word left = rotate_right(a, 2) ^ rotate_right(a, 13) ^ rotate_right(a, 22);
            const Word right = rotate_right(e, 6) ^ rotate_right(e, 11) ^ rotate_right(e, 25);
            const Word first = h + right + choose + rounds[i] + schedule[i];
            const Word second = left + majority;
            h = g;
            g = f;
            f = e;
            e = d + first;
            d = c;
            c = b;
            b = a;
            a = first + second;
        }
        const std::array<Word, 8> working{a, b, c, d, e, f, g, h};
        for (std::size_t i = 0; i < state_.size(); ++i) state_[i] += working[i];
    }

    std::array<Word, 8> state_{0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                                0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    std::array<Byte, 64> buffer_;
    std::size_t used_ = 0;
    Size length_ = 0;
};
#endif

constexpr Size chunk_size(std::string_view tag, Size payload) {
    return 1 + tag.size() + 8 + payload;
}

class FingerprintStream {
public:
    void header(std::string_view tag, Size payload) {
        byte(static_cast<Byte>(tag.size()));
        write(tag);
        little_endian(payload);
    }

    void chunk(std::string_view tag, std::string_view payload) {
        header(tag, payload.size());
        write(payload);
    }

    void byte(Byte value) {
        if (used_ == buffer_.size()) flush();
        buffer_[used_++] = value;
    }

    void little_endian(Size value) {
        if (buffer_.size() - used_ < 8) flush();
        for (unsigned i = 0; i < 8; ++i) buffer_[used_ + i] = static_cast<Byte>(value >> (8 * i));
        used_ += 8;
    }

    void number(double value, bool preserve_signed_zero = false) {
        Size bits = 0;
        if (preserve_signed_zero || value != 0) std::memcpy(&bits, &value, sizeof(bits));
        little_endian(bits);
    }

    void text(std::string_view value) { write(value); }

    HexDigest finish() {
        flush();
        return hash_.finish();
    }

private:
    void flush() {
        hash_.update(buffer_.data(), used_);
        used_ = 0;
    }

    void write(std::string_view text) {
        if (text.empty()) return;
        if (text.size() > buffer_.size() - used_) flush();
        if (text.size() >= buffer_.size()) {
            hash_.update(text.data(), text.size());
        } else {
            std::memcpy(buffer_.data() + used_, text.data(), text.size());
            used_ += text.size();
        }
    }

    Sha256 hash_;
    std::array<Byte, 16 * 1024> buffer_;
    std::size_t used_ = 0;
};

bool numeric_storage(SEXP values) {
    return TYPEOF(values) == REALSXP || TYPEOF(values) == INTSXP;
}

class NumericValues {
public:
    explicit NumericValues(SEXP values) :
        real_(TYPEOF(values) == REALSXP ? REAL(values) : nullptr),
        integer_(TYPEOF(values) == INTSXP ? INTEGER(values) : nullptr) {}

    double operator[](R_xlen_t index) const {
        if (real_) return real_[index];
        if (!integer_) return 0.0;
        const int value = integer_[index];
        return value == NA_INTEGER ? NA_REAL : static_cast<double>(value);
    }

private:
    const double* real_;
    const int* integer_;
};

struct RowSelection {
    R_xlen_t size;
    const int* rows;

    R_xlen_t operator[](R_xlen_t index) const {
        if (!rows) return index;
        return rows[index] > 0 ? rows[index] - 1 : -1;
    }
};

std::string_view scalar_tag(SEXP argument) {
    if (TYPEOF(argument) != STRSXP || XLENGTH(argument) != 1 || STRING_ELT(argument, 0) == NA_STRING) {
        Rf_error("Sample fingerprint tags must be strings.");
    }
    return CHAR(STRING_ELT(argument, 0));
}

void check_interrupt(R_xlen_t index) {
    if ((index & 0xffff) == 0) R_CheckUserInterrupt();
}

Size row_bytes(const int* positions, const RowSelection& selected, R_xlen_t n) {
    Size bytes = 0;
    for (R_xlen_t i = 0; i < selected.size; ++i) {
        check_interrupt(i);
        if (positions[i] <= 0 || selected[i] < 0 || selected[i] >= n) {
            Rf_error("Invalid sample fingerprint row selection.");
        }
        int label = positions[i] - 1;
        Size digits = 1;
        while (label >= 10) {
            label /= 10;
            ++digits;
        }
        bytes += chunk_size("int", digits);
    }
    return bytes;
}

Size column_bytes(SEXP column_names) {
    Size bytes = 0;
    for (R_xlen_t j = 0; j < XLENGTH(column_names); ++j) {
        void* temporary = vmaxget();
        const auto size = chunk_size("name", std::strlen(Rf_translateCharUTF8(STRING_ELT(column_names, j))));
        vmaxset(temporary);
        if (size > std::numeric_limits<Size>::max() - bytes) {
            Rf_error("Sample fingerprint column names exceed the length range.");
        }
        bytes += size;
    }
    return bytes;
}

void write_rows(FingerprintStream& stream, const int* positions, R_xlen_t n, Size bytes) {
    stream.header("rows", bytes);
    for (R_xlen_t i = 0; i < n; ++i) {
        check_interrupt(i);
        char decimal[std::numeric_limits<int>::digits10 + 1];
        const auto converted = std::to_chars(decimal, decimal + sizeof(decimal), positions[i] - 1);
        stream.chunk("int", {decimal, static_cast<std::size_t>(converted.ptr - decimal)});
    }
}

void write_columns(FingerprintStream& stream, SEXP column_names, Size bytes) {
    stream.header("columns", bytes);
    for (R_xlen_t j = 0; j < XLENGTH(column_names); ++j) {
        void* temporary = vmaxget();
        stream.chunk("name", Rf_translateCharUTF8(STRING_ELT(column_names, j)));
        vmaxset(temporary);
    }
}

Size vector_bytes(R_xlen_t n) {
    return 1 + 8 + sizeof(double) * static_cast<Size>(n);
}

bool write_vector(FingerprintStream& stream, std::string_view tag, const NumericValues& values,
                  const RowSelection& rows) {
    stream.header(tag, vector_bytes(rows.size));
    stream.byte(1);
    stream.little_endian(static_cast<Size>(rows.size));
    for (R_xlen_t i = 0; i < rows.size; ++i) {
        check_interrupt(i);
        const double value = values[rows[i]];
        if (!std::isfinite(value)) return false;
        stream.number(value);
    }
    return true;
}

bool write_design(FingerprintStream& stream, const NumericValues& x, R_xlen_t n,
                  const int* columns, R_xlen_t width, const RowSelection& rows) {
    stream.header("design", vector_bytes(rows.size * width) + 8);
    stream.byte(2);
    stream.little_endian(static_cast<Size>(rows.size));
    stream.little_endian(static_cast<Size>(width));
    for (R_xlen_t i = 0; i < rows.size; ++i) {
        check_interrupt(i);
        const auto row = rows[i];
        for (R_xlen_t j = 0; j < width; ++j) {
            const int column = columns[j];
            const double value = column == 0 ? 1.0 : x[(column - 1) * n + row];
            if (!std::isfinite(value)) return false;
            stream.number(value);
        }
    }
    return true;
}

}  // namespace

SEXP pm_sample_fingerprint(SEXP y, SEXP x, SEXP positions, SEXP rows, SEXP columns,
                           SEXP names, SEXP weights, SEXP weight_kind, SEXP offset,
                           SEXP offset_kind, SEXP exhausted) {
    const auto shape = Rf_getAttrib(x, R_DimSymbol);
    const R_xlen_t n = XLENGTH(y);
    if (!numeric_storage(y) || !numeric_storage(x) || TYPEOF(shape) != INTSXP ||
        XLENGTH(shape) != 2 || INTEGER(shape)[0] != n || INTEGER(shape)[1] < 0 ||
        static_cast<Size>(XLENGTH(x)) != static_cast<Size>(n) * static_cast<Size>(INTEGER(shape)[1]) ||
        TYPEOF(positions) != INTSXP || (rows != R_NilValue && TYPEOF(rows) != INTSXP) ||
        TYPEOF(columns) != INTSXP || TYPEOF(names) != STRSXP || XLENGTH(columns) != XLENGTH(names) ||
        TYPEOF(exhausted) != LGLSXP || XLENGTH(exhausted) != 1 || LOGICAL(exhausted)[0] == NA_LOGICAL) {
        Rf_error("Invalid sample fingerprint inputs.");
    }
    const auto weight_tag = scalar_tag(weight_kind);
    const auto index_tag = scalar_tag(offset_kind);
    if ((weights == R_NilValue && weight_tag != "none") ||
        (weights != R_NilValue && (!numeric_storage(weights) || XLENGTH(weights) != n ||
         (weight_tag != "analytic" && weight_tag != "frequency" && weight_tag != "probability"))) ||
        (offset != R_NilValue && (!numeric_storage(offset) || XLENGTH(offset) != n)) ||
        (index_tag != "offset" && index_tag != "exposure")) {
        Rf_error("Invalid sample fingerprint weights or offset.");
    }
    const RowSelection selected{rows == R_NilValue ? n : XLENGTH(rows),
                                rows == R_NilValue ? nullptr : INTEGER(rows)};
    const auto width = XLENGTH(columns);
    if (XLENGTH(positions) != selected.size ||
        (width > 0 && selected.size > R_XLEN_T_MAX / width)) {
        Rf_error("Invalid sample fingerprint dimensions.");
    }
    const int* column = INTEGER(columns);
    for (R_xlen_t j = 0; j < width; ++j) {
        if (column[j] < 0 || column[j] > INTEGER(shape)[1]) {
            Rf_error("Invalid sample fingerprint column selection.");
        }
        if (STRING_ELT(names, j) == NA_STRING) return Rf_ScalarString(NA_STRING);
    }

    const int* position = INTEGER(positions);
    const auto rows_size = row_bytes(position, selected, n);
    const auto names_size = column_bytes(names);
    FingerprintStream stream;
    stream.chunk("version", "powermean.sample/1");
    write_rows(stream, position, selected.size, rows_size);
    write_columns(stream, names, names_size);
    if (!write_vector(stream, "y", NumericValues(y), selected) ||
        !write_design(stream, NumericValues(x), n, column, width, selected)) {
        return Rf_ScalarString(NA_STRING);
    }
    if (weights == R_NilValue) {
        stream.chunk("weights", "none");
    } else {
        stream.header("weights", chunk_size("kind", weight_tag.size()) +
                                 chunk_size("values", vector_bytes(selected.size)));
        stream.chunk("kind", weight_tag);
        if (!write_vector(stream, "values", NumericValues(weights), selected)) {
            return Rf_ScalarString(NA_STRING);
        }
    }
    if (!write_vector(stream, index_tag, NumericValues(offset), selected)) return Rf_ScalarString(NA_STRING);
    auto digest = stream.finish();
    if (LOGICAL(exhausted)[0]) {
        Sha256 envelope;
        constexpr char tag[] = "powermean.sample/exhausted";
        envelope.update(tag, sizeof(tag));  // Include the terminating NUL.
        envelope.update(digest.data(), digest.size() - 1);
        digest = envelope.finish();
    }
    return Rf_mkString(digest.data());
}

SEXP pm_jochmans_fingerprint(SEXP prefix, SEXP codes, SEXP y, SEXP x) {
    const R_xlen_t n = XLENGTH(y);
    const auto shape = Rf_getAttrib(x, R_DimSymbol);
    if (TYPEOF(codes) != INTSXP || XLENGTH(codes) != 2 * n ||
        TYPEOF(y) != REALSXP || TYPEOF(x) != REALSXP || TYPEOF(shape) != INTSXP ||
        XLENGTH(shape) != 2 || INTEGER(shape)[0] != n ||
        static_cast<Size>(XLENGTH(x)) != static_cast<Size>(n) * static_cast<Size>(INTEGER(shape)[1])) {
        Rf_error("Invalid Jochmans sample fingerprint inputs.");
    }
    FingerprintStream stream;
    stream.text(scalar_tag(prefix));
    const int* cell = INTEGER(codes);
    const double* response = REAL(y);
    const double* design = REAL(x);
    for (R_xlen_t i = 0; i < XLENGTH(codes); ++i) {
        check_interrupt(i);
        stream.little_endian(static_cast<Size>(cell[i]));
    }
    for (R_xlen_t i = 0; i < n; ++i) {
        check_interrupt(i);
        stream.number(response[i], true);
    }
    const int width = INTEGER(shape)[1];
    for (R_xlen_t i = 0; i < n; ++i) {
        check_interrupt(i);
        for (int j = 0; j < width; ++j) stream.number(design[i + n * j], true);
    }
    return Rf_mkString(stream.finish().data());
}
