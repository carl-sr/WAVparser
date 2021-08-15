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
}

int WAV_t::write_data()
{
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
}

std::vector<uint8_t> &WAV_t::get_fmt()
{
}

std::vector<uint8_t> &WAV_t::get_data()
{
}

RIFF_t &WAV_t::get_riff()
{
}

int WAV_t::sample_size()
{
}

int WAV_t::get_sample(int i, int channel)
{
}

uint32_t WAV_t::calculate_byte_rate()
{
}

uint16_t WAV_t::calculate_block_align()
{
}

void WAV_t::clear_data()
{
}

void WAV_t::set_file_path(std::string new_file_path)
{
}

int WAV_t::write()
{
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