#include <string>
#include <cstdint>
#include <vector>
#include <limits>
#include <iostream>
#include <unordered_map>
#include <algorithm>

#include "RIFFparser.h"

#pragma once

// audio format codes - https://sites.google.com/site/musicgapi/technical-documents/wav-file-format
// only pcm and float are addressed in code, so far...
#define FORMAT_NONE 0x0000            // unknown
#define FORMAT_PCM 0x0001             // pcm/uncompressed
#define FORMAT_FLOAT 0x0003           // float encoded samples from Audacity have an audio format code of 3, don't actually know what this is but it works
#define FORMAT_MS_ADPCM 0x0002        // Microsoft ADPCM
#define FORMAT_ITU_G711_a_law 0x0006  // ITU G .711 a - law
#define FORMAT_ITU_G711_Au_law 0x0007 // ITU G .711 Âµ - law
#define FORMAT_IMA_ADPCM 0x0011       // IMA ADPCM
#define FORMAT_ITU_G723_ADPCM 0x0016  // ITU G .723 ADPCM(Yamaha)
#define FORMAT_GSM 0x0031             // GSM 6.10
#define FORMAT_ITU_G721_ADPCM 0x0040  // ITU G .721 ADPCM
#define FORMAT_MPEG 0x0050            // MPEG
#define FORMAT_EXPERIMENTAL 0xFFFF    // Experimental

#pragma pack(push, 1)
/**
 * WAV file header struct.
 */
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

#pragma pack(pop)

enum class WAV_encoding
{
    signed_16_PCM,
    signed_24_PCM,
    signed_32_PCM,
    unsigned_8_PCM,
    float_32,
    float_64,
    none
};

/**
 * Class for storing and manipulating WAV file data.
 */
template <class Sample_Type>
class WAV_t
{
public:
    struct cue_point
    {
        uint32_t sample_offset{0};
        std::string label;
    };

    struct WAV_Encoding_String
    {
        const char *encoding_str;
        WAV_encoding encoding_type;
    };

private:
#pragma pack(push, 1)
    struct cue_point_t
    {
        uint32_t identifier;
        uint32_t position;
        uint32_t data_chunk_id;
        uint32_t chunk_start;
        uint32_t block_start;
        uint32_t sample_start;
    };
#pragma pack(pop)

    // quick access
    RIFF_chunk_data_t *m_data(RIFF_chunk_list_t &data);
    RIFF_chunk_data_t *m_fmt(RIFF_chunk_list_t &fmt);

    // write fmt and data sections into RIFF_t object
    int write_fmt(RIFF_t &riff);
    int write_data(RIFF_t &riff);
    int write_cue(RIFF_t &riff);

    // Load raw byte data from the RIFF_t object
    void load_fmt(RIFF_t &riff);
    void load_data(RIFF_t &riff);
    void load_cue(RIFF_t &riff);

    // load from riff
    template <typename T>
    void load_sample_buffer_int(std::vector<uint8_t> &bytes);

    void load_sample_buffer_i24(std::vector<uint8_t> &bytes);

    template <typename T>
    void load_sample_buffer_float(std::vector<uint8_t> &bytes);

    // write to riff
    template <typename T>
    void write_sample_buffer_int(std::vector<uint8_t> &bytes);

    void write_sample_buffer_i24(std::vector<uint8_t> &bytes);

    template <typename T>
    void write_sample_buffer_float(std::vector<uint8_t> &bytes);

    // update header information based on samples vector
    void update_header();

    // generic mapping function
    template <typename From, typename To>
    static To value_map(From value, From from_min, From from_max, To to_min, To to_max);

    std::vector<std::vector<Sample_Type>> m_samples;
    WAV_fmt_t m_header;
    WAV_encoding m_encoding{WAV_encoding::none};

    std::vector<cue_point> m_cue_points;

public:
    WAV_t();

    WAV_t(std::string filename);

    int write(std::string filepath);

    // =========== METHODS FOR HEADER INFORMATION ===========

    uint16_t get_sample_rate() const;

    void set_sample_rate(uint16_t new_rate);

    void convert_sample_rate(uint16_t new_rate);

    int sample_size() const;

    std::vector<uint8_t> &extra_params();

    void print_header();

    uint32_t calculate_byte_rate();

    uint16_t calculate_block_align();

    WAV_encoding get_encoding() const;

    const std::vector<WAV_Encoding_String> &get_available_encodings() const;

    void set_encoding(WAV_encoding new_encoding);

    std::vector<WAV_t::cue_point> &cues();

    // =========== METHODS FOR SAMPLE DATA ===========

    uint16_t get_num_channels() const;

    void set_num_channels(int num_channels);

    void convert_to_mono();

    int get_num_samples();

    std::vector<Sample_Type> &channel(int i);

    std::vector<Sample_Type> &add_channel();

    void remove_channel(int i);

    void swap_channels(int a, int b);

    void clear_data();

    void reset_channel_lengths(Sample_Type fill = 0.0f);
};

// IMPLEMENTATIONS ==========================================================================================================================

// =========== PRIVATE METHODS ===========

// quick access to riff data chunk
template <class Sample_Type>
RIFF_chunk_data_t *WAV_t<Sample_Type>::m_data(RIFF_chunk_list_t &riff)
{
    return dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("data"));
}

// quick access to riff fmt chunk
template <class Sample_Type>
RIFF_chunk_data_t *WAV_t<Sample_Type>::m_fmt(RIFF_chunk_list_t &riff)
{
    return dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("fmt "));
}

template <class Sample_Type>
int WAV_t<Sample_Type>::write_fmt(RIFF_t &riff)
{
    int bytes_written{0};
    calculate_byte_rate();
    calculate_block_align();

    // directly write all bytes to a vector
    const uint8_t *fmt_bytes = reinterpret_cast<const uint8_t *>(&m_header);
    std::vector<uint8_t> bytes(16, 0);
    memcpy(&bytes.front(), fmt_bytes, 16);
    bytes_written += 16;

    // add extra params if they exist

    m_header.extra_params_size = m_header.extra_params.size();
    if (m_header.extra_params_size > 0)
    {
        bytes.push_back(reinterpret_cast<const uint8_t *>(&m_header.extra_params_size)[0]);
        bytes.push_back(reinterpret_cast<const uint8_t *>(&m_header.extra_params_size)[1]);
        bytes.insert(bytes.end(), m_header.extra_params.begin(), m_header.extra_params.end());
        bytes_written += 2 + m_header.extra_params_size;
    }

    // insert 'fmt ' chunk into riff
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("fmt ", bytes));

    return bytes_written;
}

template <class Sample_Type>
int WAV_t<Sample_Type>::write_data(RIFF_t &riff)
{
    if (m_samples.size() == 0)
        return 0;

    reset_channel_lengths();

    int bytes_written{0};

    // reserve right amount of space for storing the data
    std::vector<uint8_t> bytes;
    bytes.reserve(m_header.num_channels * m_samples[0].size());

    switch (m_encoding)
    {
    case WAV_encoding::signed_16_PCM:
        write_sample_buffer_int<int16_t>(bytes);
        break;
    case WAV_encoding::signed_24_PCM:
        write_sample_buffer_i24(bytes);
        break;
    case WAV_encoding::signed_32_PCM:
        write_sample_buffer_int<int32_t>(bytes);
        break;
    case WAV_encoding::unsigned_8_PCM:
        write_sample_buffer_int<uint8_t>(bytes);
        break;
    case WAV_encoding::float_32:
        write_sample_buffer_float<float>(bytes);
        break;
    case WAV_encoding::float_64:
        write_sample_buffer_float<double>(bytes);
        break;
    default:
        throw std::runtime_error("Unsupported audio encoding format. (" + std::to_string((int)m_encoding) + ")");
    }

    // insert data
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("data", bytes));

    return bytes_written;
}

template <class Sample_Type>
int WAV_t<Sample_Type>::write_cue(RIFF_t &riff)
{
    RIFF_chunk_list_t label_chunk("adtl"); // holds an id and an ascii string for each cue point
    std::vector<uint8_t> cue_chunk_data;   // holds cue_point_t data for each cue point

    uint32_t starting_id{1};

    // number of cue points as first four bytes
    uint32_t size = m_cue_points.size();
    for (int i = 0; i < sizeof(size); i++)
        cue_chunk_data.push_back(reinterpret_cast<uint8_t *>(&size)[i]);

    for (const auto &i : m_cue_points)
    {
        cue_point_t point = cue_point_t{
            starting_id,
            i.sample_offset,
            0x61746164, // 'data'
            0,
            0,
            i.sample_offset};

        // insert to data vector
        for (int i = 0; i < sizeof(cue_point_t); i++)
            cue_chunk_data.push_back(reinterpret_cast<uint8_t *>(&point)[i]);

        // create label data - reserve with size for string + null terminator + labl id
        std::vector<uint8_t> label_data;
        label_data.reserve(i.label.size() + 5);

        // label id to bytes
        for (int i = 0; i < 4; i++)
            label_data.push_back(reinterpret_cast<uint8_t *>(&starting_id)[i]);

        // label string to bytes
        std::for_each(i.label.begin(), i.label.end(), [&label_data](auto c)
                      { label_data.push_back(static_cast<uint8_t>(c)); });
        label_data.push_back(0); // terminator

        label_chunk.get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("labl", label_data)); // insert new labl chunks

        if (++starting_id == 0) // begins at one and increases, anything below that indicates an overflow
            throw std::runtime_error("Integer overflow on cue point id, too many cue points?");
    }

    // insert cue and label chunks
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("cue ", cue_chunk_data));
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_list_t>(std::move(label_chunk)));

    return 0;
}

template <class Sample_Type>
void WAV_t<Sample_Type>::load_fmt(RIFF_t &riff)
{
    // WAV file should have a 'fmt ' chunk - throw if fmt is not present
    RIFF_chunk_data_t *fmt = m_fmt(riff.get_root_chunk());
    if (!fmt)
        throw std::runtime_error("File does not have a valid 'fmt ' chunk.");

    memcpy(reinterpret_cast<uint8_t *>(&m_header), &fmt->get_data().front(), fmt->size());

    // load extra params if present
    if (fmt->get_data().size() > 16)
    {
        // grab size of extra params - magic value 16: end of normal header data
        mempcpy(reinterpret_cast<uint8_t *>(&m_header.extra_params_size), &fmt->get_data()[16], 2);

        // grab extra params - magic value 18: end of normal header data and size of extra params
        m_header.extra_params.reserve(m_header.extra_params_size);
        m_header.extra_params.insert(m_header.extra_params.end(), fmt->get_data().begin() + 18, fmt->get_data().end());
    }
}

template <class Sample_Type>
void WAV_t<Sample_Type>::load_data(RIFF_t &riff)
{
    // WAV file should have a 'data' chunk - throw if data is not present
    RIFF_chunk_data_t *data = m_data(riff.get_root_chunk());
    if (data == nullptr)
        throw std::runtime_error("File does not have a valid 'data' chunk.");

    std::vector<uint8_t> &d = data->get_data();

    // determine size for sample vector
    // create channels
    m_samples.insert(m_samples.begin(), m_header.num_channels, std::vector<Sample_Type>());
    // reserve space for samples
    int samples_per_channel{data->size() / m_header.num_channels};
    for (auto i : m_samples)
        i.reserve(samples_per_channel);

    // choose encoding type - none is the default
    m_encoding = WAV_encoding::none;
    if (m_header.audio_format == FORMAT_PCM) // PCM
    {
        switch (m_header.bits_per_sample)
        {
        case 8:
            m_encoding = WAV_encoding::unsigned_8_PCM;
            break;
        case 16:
            m_encoding = WAV_encoding::signed_16_PCM;
            break;
        case 24:
            m_encoding = WAV_encoding::signed_24_PCM;
            break;
        case 32:
            m_encoding = WAV_encoding::signed_32_PCM;
            break;
        }
    }
    else if (m_header.audio_format == FORMAT_FLOAT) // float
    {
        switch (m_header.bits_per_sample)
        {
        case 32:
            m_encoding = WAV_encoding::float_32;
            break;
        case 64:
            m_encoding = WAV_encoding::float_64;
            break;
        }
    }

    // choose read function based on encoding type
    switch (m_encoding)
    {
    case WAV_encoding::signed_16_PCM:
        load_sample_buffer_int<int16_t>(d);
        break;
    case WAV_encoding::signed_24_PCM:
        load_sample_buffer_i24(d);
        break;
    case WAV_encoding::signed_32_PCM:
        load_sample_buffer_int<int32_t>(d);
        break;
    case WAV_encoding::unsigned_8_PCM:
        load_sample_buffer_int<uint8_t>(d);
        break;
    case WAV_encoding::float_32:
        load_sample_buffer_float<float>(d);
        break;
    case WAV_encoding::float_64:
        load_sample_buffer_float<double>(d);
        break;
    default:
        throw std::runtime_error("Unsupported audio encoding format. (" + std::to_string((int)m_encoding) + ")");
    }
}

template <class Sample_Type>
void WAV_t<Sample_Type>::load_cue(RIFF_t &riff)
{
    std::unordered_map<uint32_t, std::string> labl_identifiers;

    // find labels - these are optional?
    for (auto &i : riff.get_root_chunk().get_subchunks())
    {
        // find and search list chunks
        RIFF_chunk_list_t *list = dynamic_cast<RIFF_chunk_list_t *>(i.get());
        if (list != nullptr)
        {
            // looking for 'associated data list'
            if (strcmp(list->get_form_type(), "adtl") != 0)
                continue;

            // this is the correct list chunk - find all of the 'labl' or 'note' chunks - they contain the label strings
            const std::vector<std::unique_ptr<RIFF_chunk_t>> &adtl_chunks = list->get_subchunks();
            for (auto &j : adtl_chunks)
            {
                RIFF_chunk_data_t *adtl = dynamic_cast<RIFF_chunk_data_t *>(j.get());
                if (adtl != nullptr && (strcmp(adtl->get_identifier(), "labl") == 0 || strcmp(adtl->get_identifier(), "note") == 0))
                {
                    // this is an 'adtl' or 'note' chunk
                    uint32_t adtl_id;
                    memcpy(&adtl_id, &adtl->get_data().front(), sizeof(adtl_id));
                    char *identifier = reinterpret_cast<char *>(&adtl->get_data().front()) + 4;
                    labl_identifiers[adtl_id] = std::string(identifier, strlen(identifier));
                }
            }
        }
    }

    // grab cue
    RIFF_chunk_data_t *cue_chunk = dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("cue "));
    if (cue_chunk != nullptr)
    {
        std::vector<uint8_t> cue_v = cue_chunk->get_data();

        // copy length
        uint32_t cue_chunk_length{0};
        memcpy(&cue_chunk_length, &cue_v.front(), sizeof(cue_chunk_length));
        cue_v.erase(cue_v.begin(), cue_v.begin() + sizeof(cue_chunk_length)); // delete this data

        // grab each cue chunk
        for (int i = 0; i < cue_chunk_length; i++)
        {
            cue_point_t cue_point;
            memcpy(&cue_point, &cue_v.front(), sizeof(cue_point));
            cue_v.erase(cue_v.begin(), cue_v.begin() + sizeof(cue_point));

            WAV_t<Sample_Type>::cue_point point{
                cue_point.sample_start,
                labl_identifiers.find(cue_point.identifier) == labl_identifiers.end() ? std::to_string(cue_point.identifier) : labl_identifiers[cue_point.identifier]};

            m_cue_points.push_back(point);
        }
    }
}

template <class Sample_Type>
template <class T>
void WAV_t<Sample_Type>::load_sample_buffer_int(std::vector<uint8_t> &bytes)
{
    T *buffer = reinterpret_cast<T *>(&bytes.front());
    int channel_counter = 0;
    int total_samples = bytes.size() / sizeof(T);

    // assign each sample to a channel
    for (int i = 0; i < total_samples; i++)
    {
        Sample_Type new_value = value_map<T, Sample_Type>(buffer[i], std::numeric_limits<T>::min(), std::numeric_limits<T>::max(), -1.0, 1.0);
        m_samples[channel_counter++ % m_header.num_channels].push_back(new_value);
    }
}

template <class Sample_Type>
void WAV_t<Sample_Type>::load_sample_buffer_i24(std::vector<uint8_t> &bytes)
{
    uint8_t *buffer = reinterpret_cast<uint8_t *>(&bytes.front());
    int channel_counter = 0;
    int total_samples = bytes.size() / 3; // size of 24 bit integer = 3

    // assign each sample to a channel
    for (int i = 0; i < total_samples; i += 3)
    {
        int32_t smp = buffer[i] + (buffer[i + 1] << 8) + (buffer[i + 2] << 16);

        // 24i max:  0x7fffff
        // 24i min: -0x800000
        Sample_Type new_value = value_map<int32_t, Sample_Type>(smp, -0x800000, 0x7fffff, -1.0, 1.0);
        m_samples[channel_counter++ % m_header.num_channels].push_back(new_value);
    }
}

template <class Sample_Type>
template <class T>
void WAV_t<Sample_Type>::load_sample_buffer_float(std::vector<uint8_t> &bytes)
{
    T *buffer = reinterpret_cast<T *>(&bytes.front());
    int channel_counter = 0;
    int total_samples = bytes.size() / sizeof(T);

    // assign each sample to a channel
    for (int i = 0; i < total_samples; i++)
        m_samples[channel_counter++ % m_header.num_channels].push_back(static_cast<Sample_Type>(buffer[i]));
}

template <class Sample_Type>
template <class From, class To>
To WAV_t<Sample_Type>::value_map(From value, From from_min, From from_max, To to_min, To to_max)
{
    Sample_Type d_value = static_cast<Sample_Type>(value);
    Sample_Type d_from_min = static_cast<Sample_Type>(from_min);
    Sample_Type d_from_max = static_cast<Sample_Type>(from_max);
    Sample_Type d_to_min = static_cast<Sample_Type>(to_min);
    Sample_Type d_to_max = static_cast<Sample_Type>(to_max);

    if (d_value >= d_from_max)
        return to_max;

    if (d_value <= d_from_min)
        return to_min;

    return static_cast<To>((d_value - d_from_min) * (d_to_max - d_to_min) / (d_from_max - d_from_min) + d_to_min);
}

template <class Sample_Type>
template <class T>
void WAV_t<Sample_Type>::write_sample_buffer_int(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : m_samples)
        total_samples += i.size();
    bytes.reserve(total_samples * sizeof(T));

    for (int i = 0; i < total_samples; i++)
    {
        T smp = value_map<Sample_Type, T>((m_samples[i % m_header.num_channels][i / m_header.num_channels]), -1.0, 1.0, std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

        // push individual bytes from smp into bytes
        for (int j = 0; j < sizeof(T); j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

template <class Sample_Type>
void WAV_t<Sample_Type>::write_sample_buffer_i24(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : m_samples)
        total_samples += i.size();
    bytes.reserve(total_samples * 3);

    for (int i = 0; i < total_samples; i++)
    {
        int32_t smp = value_map<Sample_Type, int32_t>(m_samples[i % m_header.num_channels][i / m_header.num_channels], -1.0, 1.0, -0x800000, 0x7fffff);

        // push individual bytes from smp into bytes
        for (int j = 0; j < 3; j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

template <class Sample_Type>
template <class T>
void WAV_t<Sample_Type>::write_sample_buffer_float(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : m_samples)
        total_samples += i.size();
    bytes.reserve(total_samples * sizeof(T));

    for (int i = 0; i < total_samples; i++)
    {
        T smp = static_cast<T>(m_samples[i % m_header.num_channels][i / m_header.num_channels]);

        // push individual bytes from smp into bytes
        for (int j = 0; j < sizeof(T); j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

template <class Sample_Type>
void WAV_t<Sample_Type>::update_header()
{
    m_header.num_channels = m_samples.size();
    switch (m_encoding)
    {
    case WAV_encoding::signed_16_PCM:
        m_header.audio_format = FORMAT_PCM;
        m_header.bits_per_sample = 16;
        break;
    case WAV_encoding::signed_24_PCM:
        m_header.audio_format = FORMAT_PCM;
        m_header.bits_per_sample = 24;
        break;
    case WAV_encoding::signed_32_PCM:
        m_header.audio_format = FORMAT_PCM;
        m_header.bits_per_sample = 32;
        break;
    case WAV_encoding::unsigned_8_PCM:
        m_header.audio_format = FORMAT_PCM;
        m_header.bits_per_sample = 8;
        break;
    case WAV_encoding::float_32:
        m_header.audio_format = FORMAT_FLOAT;
        m_header.bits_per_sample = 32;
        break;
    case WAV_encoding::float_64:
        m_header.audio_format = FORMAT_FLOAT;
        m_header.bits_per_sample = 64;
        break;
    }
    calculate_block_align();
    calculate_byte_rate();
}

// =========== PUBLIC METHODS ===========

/**
 * Construct an empty WAV file. Contains no samples and 
 * default header data.
 */
template <class Sample_Type>
WAV_t<Sample_Type>::WAV_t() : m_samples(2, std::vector<Sample_Type>())
{
    m_encoding = WAV_encoding::signed_32_PCM;
}

/**
 * Construct a WAV_t object from a WAV file.
 * @param filename The file to parse.
 */
template <class Sample_Type>
WAV_t<Sample_Type>::WAV_t(std::string filename)
{
    RIFF_t riff(filename);

    if (strcmp(riff.get_root_chunk().get_form_type(), "WAVE") != 0)
        throw std::runtime_error("File is not a valid WAVE file.");

    load_fmt(riff);
    load_data(riff);
    load_cue(riff);
}

/**
 * Write WAV_t data to the disk at the specified filepath. 
 * @return Number of bytes written
 */
template <class Sample_Type>
int WAV_t<Sample_Type>::write(std::string filepath)
{
    update_header();
    RIFF_t riff;
    riff.get_root_chunk().set_form_type("WAVE");
    write_fmt(riff);
    write_data(riff);
    write_cue(riff);
    riff.set_filepath(filepath);

    return riff.write();
}

// =========== METHODS FOR HEADER INFORMATION ===========

/**
 * Get the current sample rate of the WAV data.
 */
template <class Sample_Type>
uint16_t WAV_t<Sample_Type>::get_sample_rate() const
{
    return m_header.sample_rate;
}

/**
 * Set a new sample rate for the WAV data.
 * @param new_rate The new sample rate.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::set_sample_rate(uint16_t new_rate)
{
    m_header.sample_rate = new_rate;
}

/**
 * Convert sample rate without affecting audio (too  much).
 * @param new_rate The new sample rate.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::convert_sample_rate(uint16_t new_rate)
{
    // TODO
}

/**
 * Get the size in bytes of a sample when written.
 */
template <class Sample_Type>
int WAV_t<Sample_Type>::sample_size() const
{
    switch (m_encoding)
    {
    case WAV_encoding::signed_16_PCM:
        return 2;
    case WAV_encoding::signed_24_PCM:
        return 3;
    case WAV_encoding::signed_32_PCM:
        return 4;
    case WAV_encoding::unsigned_8_PCM:
        return 1;
    case WAV_encoding::float_32:
        return 4;
    case WAV_encoding::float_64:
        return 8;
    default:
        return 0;
    };
}

/**
 * Get the vector containing byte data for extra parameters. 
 * Extra parameters are not parsed by the object.
 */
template <class Sample_Type>
std::vector<uint8_t> &WAV_t<Sample_Type>::extra_params()
{
    return m_header.extra_params;
}

/**
 * Quickly print header information
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::print_header()
{
    printf("audio format: %d\n", m_header.audio_format);
    printf("num channels: %d\n", m_header.num_channels);
    printf("sample rate %d\n", m_header.sample_rate);
    printf("byte rate: %d\n", m_header.byte_rate);
    printf("block align %d\n", m_header.block_align);
    printf("bits per sample %d\n", m_header.bits_per_sample);

    if (m_header.extra_params_size > 0)
    {
        printf("extra params size: %d, extra params:\n", m_header.extra_params_size);
        for (auto i : m_header.extra_params)
            printf(" %02x", i);
        putchar('\n');
    }
}

/**
 * Helper function to set header byte rate. Generally only used internally but can be useful to do 
 * correct calculations when changing header information manually.
 * @return The new calculated byte rate.
 */
template <class Sample_Type>
uint32_t WAV_t<Sample_Type>::calculate_byte_rate()
{
    m_header.byte_rate = m_header.sample_rate * m_header.num_channels * sample_size();
    return m_header.byte_rate;
}

/**
 * Helper function to set header block align. Generally only used internally but can be useful to do 
 * correct calculations when changing header information manually.
 * @return The new calculated block align.
 */
template <class Sample_Type>
uint16_t WAV_t<Sample_Type>::calculate_block_align()
{
    m_header.block_align = m_header.num_channels * sample_size();
    return m_header.block_align;
}

/**
 * Get the current encoding type of the WAV data.
 */
template <class Sample_Type>
WAV_encoding WAV_t<Sample_Type>::get_encoding() const
{
    return m_encoding;
}

/**
 * Get a list of the operational WAV encodings.
 */
template <class Sample_Type>
const std::vector<typename WAV_t<Sample_Type>::WAV_Encoding_String> &WAV_t<Sample_Type>::get_available_encodings() const
{
    static const std::vector<WAV_t<Sample_Type>::WAV_Encoding_String> available_encodings = {
        {"16-bit PCM", WAV_encoding::signed_16_PCM},
        {"24-bit PCM", WAV_encoding::signed_24_PCM},
        {"32-bit PCM", WAV_encoding::signed_32_PCM},
        {"8-bit PCM", WAV_encoding::unsigned_8_PCM},
        {"32-bit Floating Point", WAV_encoding::float_32},
        {"64-bit Floating Point", WAV_encoding::float_64},
    };

    return available_encodings;
}

/**
 * Set a new encoding type for the WAV data
 * @param new_encoding The encoding type to use during the next write.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::set_encoding(WAV_encoding new_encoding)
{
    m_encoding = new_encoding;
    update_header();
}

/**
 * Get the list of cue points.
 * @return Reference to list of cue points.
 */
template <class Sample_Type>
std::vector<typename WAV_t<Sample_Type>::cue_point> &WAV_t<Sample_Type>::cues()
{
    return m_cue_points;
}

// =========== METHODS FOR SAMPLE DATA ===========

/**
 * Get the current number of audio channels.
 */
template <class Sample_Type>
uint16_t WAV_t<Sample_Type>::get_num_channels() const
{
    return m_samples.size();
}

/**
 * Set the number of channels
 * If the new number of channels is smaller than the current number, channels will be deleted.
 * If it is greater than the current number, new channels initialized to zero will be created.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::set_num_channels(int num_channels)
{
    if (get_num_channels() == num_channels)
        return;

    while (get_num_channels() != num_channels)
    {
        if (get_num_channels() > num_channels)
            m_samples.pop_back();
        else
            m_samples.emplace_back(get_num_samples(), 0.0f);
    }
}

/**
 * Convert entire file to mono. Averages all channel samples into channel 0
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::convert_to_mono()
{
    reset_channel_lengths();
    for (int i = 0; i < get_num_samples(); i++)
    {
        Sample_Type sum{0.0f};
        std::for_each(m_samples.begin(), m_samples.end(), [i, &sum](auto a)
                      { sum += a[i]; });
        channel(0)[i] = sum / get_num_channels();
    }
}

/**
 * Get the current number of samples in the sample buffer
 */
template <class Sample_Type>
int WAV_t<Sample_Type>::get_num_samples()
{
    if (get_num_channels() == 0)
        return 0;
    return channel(0).size();
}

/**
 * Get channel data.
 * @param i The channel number to grab
 */
template <class Sample_Type>
std::vector<Sample_Type> &WAV_t<Sample_Type>::channel(int i)
{
    return m_samples[i];
}

/**
 * Add an audio channel.
 * @returns Reference to the newly created audio channel
 */
template <class Sample_Type>
std::vector<Sample_Type> &WAV_t<Sample_Type>::add_channel()
{
    m_samples.emplace_back(get_num_samples(), 0.0f);
    m_header.num_channels++;
    return m_samples.back();
}

/**
 * Remove an audio channel.
 * @param i The index of the channel to remove
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::remove_channel(int i)
{
    m_samples.erase(m_samples.begin() + i);
    m_header.num_channels--;
}

/**
 * Swap two audio channels.
 * @param a First channel to swap.
 * @param b Second channel to swap.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::swap_channels(int a, int b)
{
    std::swap(m_samples[a], m_samples[b]);
}

/**
 * Clear all sample data from the file.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::clear_data()
{
    m_samples = std::vector<std::vector<Sample_Type>>(m_header.num_channels, std::vector<Sample_Type>());
}

/**
 * Set the length of all channels to the length of the longest channel.
 * @param fill The value to append to short channels.
 */
template <class Sample_Type>
void WAV_t<Sample_Type>::reset_channel_lengths(Sample_Type fill)
{
    // find the length of the longest sample vector
    int len{0};
    for (auto i : m_samples)
        len = std::max(static_cast<int>(i.size()), len);

    // fill the end of each sample vector
    for (auto i : m_samples)
        if (i.size() < len)
            i.insert(i.end(), fill, len - i.size());
}