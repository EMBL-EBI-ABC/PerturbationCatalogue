// Validate complete SRA spots and stream barcode/biological pairs to stdout.
#include <array>
#include <charconv>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Read { uint64_t id; std::string header, sequence, quality; };

uint64_t number(const std::string& s, size_t begin, size_t end) {
    uint64_t value = 0;
    auto parsed = std::from_chars(s.data() + begin, s.data() + end, value);
    if (begin == end || parsed.ec != std::errc() || parsed.ptr != s.data() + end || !value)
        throw std::runtime_error("Invalid spot/read identifier");
    return value;
}

int main(int argc, char** argv) {
    try {
        if (argc != 5)
            throw std::runtime_error("Usage: read_router ACCESSION FASTQ WORKFLOW CHEMISTRY");
        const bool kite = std::string(argv[3]) == "kite";
        const bool v1 = std::string(argv[4]) == "10xv1";
        std::array<char, 1 << 20> input_buffer;
        std::ifstream input;
        input.rdbuf()->pubsetbuf(input_buffer.data(), input_buffer.size());
        input.open(argv[2], std::ios::binary);
        if (!input) throw std::runtime_error("Cannot open FASTQ");
        std::string prefix = "@" + std::string(argv[1]) + ".";
        std::string header, sequence, plus, quality, output;
        output.reserve(2 << 20);
        std::vector<Read> reads;
        std::vector<std::pair<uint64_t, size_t>> layout;
        uint64_t spot = 0, spots = 0, records = 0;
        auto flush = [&]() {
            if (std::fwrite(output.data(), 1, output.size(), stdout) != output.size())
                throw std::runtime_error("Counter pipe write failed");
            output.clear();
        };
        auto emit = [&]() {
            std::vector<std::pair<uint64_t, size_t>> current;
            const Read *barcode = nullptr, *biological = nullptr;
            Read combined;
            for (const auto& r : reads) {
                current.emplace_back(r.id, r.sequence.size());
                if (r.sequence.size() > 40) {
                    if (biological) throw std::runtime_error("Ambiguous biological read");
                    biological = &r;
                }
            }
            if (v1) {
                const Read *cell = nullptr, *umi = nullptr;
                std::vector<const Read*> short_reads;
                for (const auto& r : reads) {
                    if (&r == biological) continue;
                    if (r.sequence.size() >= 13 && r.sequence.size() <= 16) {
                        if (cell) throw std::runtime_error("Ambiguous 10x v1 cell-barcode read");
                        cell = &r;
                    } else if (r.sequence.size() >= 8 && r.sequence.size() <= 12) {
                        short_reads.push_back(&r);
                    }
                }
                for (const auto* r : short_reads) {
                    if (r->sequence.size() == 10) {
                        if (umi) throw std::runtime_error("Ambiguous 10x v1 UMI read");
                        umi = r;
                    }
                }
                if (!umi && short_reads.size() == 1) umi = short_reads.front();
                if (!cell || !umi) {
                    throw std::runtime_error("Missing 10x v1 cell-barcode/UMI reads");
                }
                if (kite) {
                    const Read *guide = nullptr;
                    for (const auto& r : reads) {
                        if (&r == cell || &r == umi || &r == biological) continue;
                        if (guide) throw std::runtime_error("Ambiguous 10x v1 guide read");
                        guide = &r;
                    }
                    if (guide) biological = guide;
                }
                if (!biological) throw std::runtime_error("Missing biological read");
                combined = *cell;
                combined.sequence += umi->sequence;
                combined.quality += umi->quality;
                barcode = &combined;
            } else {
                for (const auto& r : reads) {
                    if (r.sequence.size() >= 20 && r.sequence.size() <= 40) {
                        if (barcode) throw std::runtime_error("Ambiguous barcode read");
                        barcode = &r;
                    }
                }
            }
            if (!barcode || !biological) {
                throw std::runtime_error("Missing barcode or biological read");
            }
            if (layout.empty()) layout = current;
            if (current != layout) throw std::runtime_error("Read layout changed within accession");
            for (const auto* r : {barcode, biological}) {
                output.append(r->header).push_back('\n');
                output.append(r->sequence).append("\n+\n");
                output.append(r->quality).push_back('\n');
            }
            if (output.size() >= (1 << 20)) flush();
            ++spots;
            reads.clear();
        };
        while (std::getline(input, header)) {
            if (!std::getline(input, sequence) || !std::getline(input, plus) ||
                !std::getline(input, quality) || input.eof() ||
                header.compare(0, prefix.size(), prefix) || plus != "+" ||
                sequence.empty() || sequence.size() != quality.size())
                throw std::runtime_error("Malformed or truncated FASTQ record");
            const auto slash = header.find('/', prefix.size());
            if (slash == std::string::npos) throw std::runtime_error("Missing read identifier");
            auto next = number(header, prefix.size(), slash);
            auto id = number(header, slash + 1, header.size());
            if (next != spot) {
                if (spot) {
                    if (next != spot + 1) throw std::runtime_error("Missing or unordered spot");
                    emit();
                } else if (next != 1) throw std::runtime_error("First spot must be 1");
                spot = next;
            }
            if (!reads.empty() && id <= reads.back().id)
                throw std::runtime_error("Duplicate or unordered read identifier");
            reads.push_back({id, header, sequence, quality});
            ++records;
        }
        if (input.bad() || !spot) throw std::runtime_error("Unreadable or empty FASTQ");
        emit();
        flush();
        if (std::fflush(stdout)) throw std::runtime_error("Counter pipe flush failed");
        std::cerr << "{\"spots\":" << spots << ",\"input_records\":" << records << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Read routing failed: " << error.what() << '\n';
        return 1;
    }
}
