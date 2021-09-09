#include "WAVparser.h"

// =========== PIRVATE METHODS ===========
RIFF_chunk_data_t *WAV_t::m_data(RIFF_chunk_list_t &riff)
{
    return reinterpret_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("data"));
}

RIFF_chunk_data_t *WAV_t::m_fmt(RIFF_chunk_list_t &riff)
{
    return reinterpret_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("fmt "));
}

int WAV_t::write_fmt(RIFF_t &riff)
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
    if (header.extra_params_size > 0)
    {
        bytes.push_back(reinterpret_cast<const uint8_t *>(&header.extra_params_size)[0]);
        bytes.push_back(reinterpret_cast<const uint8_t *>(&header.extra_params_size)[1]);
        bytes.insert(bytes.end(), header.extra_params.begin(), header.extra_params.end());
        bytes_written += 2 + header.extra_params_size;
    }

    // insert 'fmt ' chunk into riff
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("fmt ", bytes));

    return bytes_written;
}

int WAV_t::write_data(RIFF_t &riff)
{
    if (samples.size() == 0)
        return 0;

    reset_channel_lengths();

    int bytes_written{0};

    // reserve right amount of space for storing the data
    std::vector<uint8_t> bytes;
    bytes.reserve(header.num_channels * samples[0].size());

    // TODO: convert samples into bytes

    // insert 'data' chunk into riff
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("data", bytes));

    return bytes_written;
}

void WAV_t::load_fmt(RIFF_chunk_data_t &fmt)
{
    memcpy(reinterpret_cast<uint8_t *>(&header), &fmt.get_data().front(), fmt.size());

    // load extra params if present
    if (fmt.get_data().size() > 16)
    {
        // grab size of extra params - magic value 16: end of normal header data
        mempcpy(reinterpret_cast<uint8_t *>(&header.extra_params_size), &fmt.get_data()[16], 2);

        // grab extra params - magic value 18: end of normal header data and size of extra params
        header.extra_params.reserve(header.extra_params_size);
        header.extra_params.insert(header.extra_params.end(), fmt.get_data().begin() + 18, fmt.get_data().end());
    }
}

void WAV_t::load_data(RIFF_chunk_data_t &data)
{
    std::vector<uint8_t> &d = data.get_data();

    // determine size for sample vector
    // create channels
    samples.insert(samples.begin(), header.num_channels, std::vector<float>());
    // reserve space for samples
    int samples_per_channel{data.size() / header.num_channels};
    for (auto i : samples)
        i.reserve(samples_per_channel);

    // TODO: read byte data into audio channels
}

// =========== PUBLIC METHODS ===========

WAV_t::WAV_t() : samples(2, std::vector<float>()) {}

WAV_t::WAV_t(std::string filename)
{
    RIFF_t riff(filename);

    if (strcmp(riff.get_root_chunk().get_form_type(), "WAVE") != 0)
        throw std::runtime_error("File is not a valid WAVE file.");

    // WAV file should have a 'fmt ' chunk
    RIFF_chunk_data_t *fmt_chunk = dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("fmt "));
    if (!fmt_chunk)
        throw std::runtime_error("File does not have a valid 'fmt ' chunk.");

    // WAV file should have a 'data' chunk
    RIFF_chunk_data_t *data_chunk = dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("data"));
    if (!data_chunk)
        throw std::runtime_error("File does not have a valid 'data' chunk.");

    // load required chunk data
    load_fmt(*fmt_chunk);
    load_data(*data_chunk);
}

int WAV_t::write(std::string filepath)
{
    RIFF_t riff;
    riff.get_root_chunk().set_form_type("WAVE");
    write_fmt(riff);
    write_data(riff);
    riff.set_filepath(filepath);

    return riff.write();
}
// =========== METHODS FOR HEADER INFORMATION ===========

uint16_t WAV_t::sample_rate() const
{
    return header.sample_rate;
}

void WAV_t::set_sample_rate(uint16_t new_rate)
{
    header.sample_rate = new_rate;
}

int WAV_t::sample_size() const
{
    return header.bits_per_sample / 8;
}

void WAV_t::set_sample_size(int new_size)
{
    header.bits_per_sample = new_size * 8;
}

std::vector<uint8_t> &WAV_t::extra_params()
{
    return header.extra_params;
}

void WAV_t::print_header()
{
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

WAV_encoding WAV_t::get_encoding() const
{
    return encoding;
}

void WAV_t::set_encoding(WAV_encoding new_encoding)
{
    encoding = new_encoding;
}
// =========== METHODS FOR SAMPLE DATA ===========

uint16_t WAV_t::num_channels() const
{
    return samples.size();
}

std::vector<float> &WAV_t::channel(int i)
{
    return samples[i];
}
std::vector<float> &WAV_t::add_channel()
{
    samples.push_back({});
    header.num_channels++;
    return samples.back();
}
void WAV_t::remove_channel(int i)
{
    samples.erase(samples.begin() + i);
    header.num_channels--;
}

void WAV_t::swap_channels(int a, int b)
{
    std::swap(samples[a], samples[b]);
}

void WAV_t::clear_data()
{
    samples = std::vector<std::vector<float>>(header.num_channels, std::vector<float>());
}

void WAV_t::reset_channel_lengths(float fill)
{
    // find the length of the longest sample vector
    int len{0};
    for (auto i : samples)
        len = std::max(static_cast<int>(i.size()), len);

    // fill the end of each sample vector
    for (auto i : samples)
        if (i.size() < len)
            i.insert(i.end(), fill, len - i.size());
}