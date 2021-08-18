#include "WAVparser.h"

WAV_t::WAV_t()
{
    m_riff.get_root_chunk().set_form_type("WAVE");

    // add format and data chunks
    std::vector<std::unique_ptr<RIFF_chunk_t>> &chunks = m_riff.get_root_chunk().get_subchunks();
    chunks.push_back(std::make_unique<RIFF_chunk_data_t>("fmt "));
    chunks.push_back(std::make_unique<RIFF_chunk_data_t>("data"));
}

WAV_t::WAV_t(std::string filename) : m_riff(filename)
{
    if (strcmp(m_riff.get_root_chunk().get_form_type(), "WAVE") != 0)
        throw std::runtime_error("File is not a valid WAVE file.");

    if (!m_riff.exists_chunk_with_id("fmt "))
        throw std::runtime_error("File does not have a 'fmt ' chunk.");

    if (!m_riff.exists_chunk_with_id("data"))
        throw std::runtime_error("File does not have a 'data' chunk.");

    load_fmt();
    load_data();
}

RIFF_chunk_data_t *WAV_t::m_data()
{
    return reinterpret_cast<RIFF_chunk_data_t *>(m_riff.get_chunk_with_id("data"));
}

RIFF_chunk_data_t *WAV_t::m_fmt()
{
    return reinterpret_cast<RIFF_chunk_data_t *>(m_riff.get_chunk_with_id("fmt "));
}

int WAV_t::write_fmt()
{
    int bytes_written{0};
    calculate_byte_rate();
    calculate_block_align();

    // directly write all bytes to a vector
    const uint8_t *fmt_bytes = reinterpret_cast<const uint8_t *>(&header);
    std::vector<uint8_t> bytes(16, 0);
    memcpy(&bytes.front(), fmt_bytes, 16);
    bytes_written += 16;

    // add extra params if they exist

    header.extra_params_size = header.extra_params.size();
    if(header.extra_params_size > 0)
    {
        bytes.push_back(reinterpret_cast<const uint8_t *>(&header.extra_params_size)[0]);
        bytes.push_back(reinterpret_cast<const uint8_t *>(&header.extra_params_size)[1]);
        bytes.insert(bytes.end(), header.extra_params.begin(), header.extra_params.end());
        bytes_written += 2 + header.extra_params_size;
    }

    // set new fmt data
    m_fmt()->set_data(bytes);

    return bytes_written;
}

int WAV_t::write_data()
{
    int bytes_written{0};

    // determine size of each sample
    int bytes_per_sample = header.bits_per_sample / 8;

    std::vector<uint8_t> bytes;
    bytes.reserve(bytes_per_sample * samples.size());

    for(auto i : samples)
    {
        for(int b = 0; b < bytes_per_sample; b++)
        {
            // get one byte out of each sample
            // shift to get each byte out of the uint64_t and mask with 0xff
            bytes.push_back((i >> (bytes_per_sample - b - 1) * 8) & 0xff);
            bytes_written += 1;
        }
    }

    m_data()->set_data(bytes);

    return bytes_written;
}

void WAV_t::load_fmt()
{
    memcpy(reinterpret_cast<uint8_t *>(&header), &m_fmt()->get_data().front(), m_fmt()->size());
}

void WAV_t::load_data()
{
    std::vector<uint8_t> &d = m_data()->get_data();

    // determine size for sample vector
    int bytes_per_sample = header.bits_per_sample / 8;
    samples.reserve(m_data()->size() / bytes_per_sample);

    uint64_t smp{0};
    for (int i = 0; i < d.size(); i++)
    {

        if (i % bytes_per_sample == 0 && i != 0)
        {
            samples.push_back(smp);
            smp = 0;
        }

        smp += (d[i] << ((bytes_per_sample - (i % bytes_per_sample) - 1) * 8));
    }
    // last sample
    samples.push_back(smp);
}

std::vector<uint8_t> &WAV_t::get_fmt()
{
    return m_fmt()->get_data();
}

std::vector<uint8_t> &WAV_t::get_data()
{
    return m_data()->get_data();
}

RIFF_t &WAV_t::get_riff()
{
    return m_riff;
}

int WAV_t::sample_size()
{
    return header.bits_per_sample / 8;
}

uint64_t &WAV_t::get_sample(int i, int channel)
{
    if(channel > (header.num_channels - 1))
        throw std::runtime_error("Requested access to audio channel that does not exist");

    return samples[i + channel];
}

uint32_t WAV_t::calculate_byte_rate()
{
    header.byte_rate = header.sample_rate * header.num_channels * sample_size();
    return header.byte_rate;
}

uint16_t WAV_t::calculate_block_align()
{
    header.block_align = header.num_channels * sample_size();
    return header.block_align;
}

void WAV_t::clear_data()
{
    samples.clear();
    write_data();
}

void WAV_t::set_filepath(std::string new_file_path)
{
    m_riff.set_filepath(new_file_path);
}

int WAV_t::write()
{
    write_data();
    write_fmt();

    return m_riff.write();
}

void WAV_t::print_header()
{
    printf("*** %s header ***\n", m_riff.get_filepath().c_str());
    printf("audio format: %d\n", header.audio_format);
    printf("num channels: %d\n", header.num_channels);
    printf("sample rate %d\n", header.sample_rate);
    printf("byte rate: %d\n", header.byte_rate);
    printf("block align %d\n", header.block_align);
    printf("bits per sample %d\n", header.bits_per_sample);

    if (header.extra_params_size > 0)
    {
        printf("extra params size: %d, extra params:\n", header.extra_params_size);
        for (auto i : header.extra_params)
            printf(" %02x", i);
        putchar('\n');
    }
}