#include <string>
#include <cstdint>
#include <vector>

#include "RIFFparser.h"

#pragma once

#pragma pack(2)
struct WAV_fmt_t
{
    uint16_t audio_format{1};
    uint16_t num_channels{2};
    uint32_t sample_rate{44100};
    uint32_t byte_rate{176400};
    uint16_t block_align{4};
    uint16_t bits_per_sample{16};
    uint16_t extra_params_size{0};
    std::vector<uint8_t> extra_params; // generally don't exist
};

class WAV_t
{
private:
    RIFF_t m_riff;

    // quick access
    RIFF_chunk_data_t *m_data();
    RIFF_chunk_data_t *m_fmt();

    // write fmt and data sections into m_riff
    int write_fmt();
    int write_data();

public:
    WAV_fmt_t header;
    std::vector<uint64_t> samples;

    WAV_t();
    WAV_t(std::string filename);

    void load_fmt();
    void load_data();

    std::vector<uint8_t> &get_fmt();
    std::vector<uint8_t> &get_data();
    RIFF_t &get_riff();

    int sample_size();
    uint64_t &get_sample(int i, int channel = 0);

    uint32_t calculate_byte_rate();
    uint16_t calculate_block_align();

    void clear_data();

    void set_filepath(std::string new_file_path);

    int write();

    void print_header();
};