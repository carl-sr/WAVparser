#include "WAVparser.h"

// =========== PIRVATE METHODS ===========
RIFF_chunk_data_t *WAV_t::m_data(RIFF_chunk_list_t &riff)
{
    return dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("data"));
}

RIFF_chunk_data_t *WAV_t::m_fmt(RIFF_chunk_list_t &riff)
{
    return dynamic_cast<RIFF_chunk_data_t *>(riff.get_chunk_with_id("fmt "));
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

    switch (encoding)
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
        throw std::runtime_error("Unsupported audio encoding format. (" + std::to_string((int)encoding) + ")");
    }

    // insert data
    riff.get_root_chunk().get_subchunks().push_back(std::make_unique<RIFF_chunk_data_t>("data", bytes));

    return bytes_written;
}

int WAV_t::write_cue(RIFF_t &riff)
{
    RIFF_chunk_list_t label_chunk("adtl"); // holds an id and an ascii string for each cue point
    std::vector<uint8_t> cue_chunk_data;   // holds cue_point_t data for each cue point

    uint32_t starting_id{1};

    // number of cue points as first four bytes
    uint32_t size = cue_points.size();
    for (int i = 0; i < sizeof(size); i++)
        cue_chunk_data.push_back(reinterpret_cast<uint8_t *>(&size)[i]);

    for (const auto &i : cue_points)
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

void WAV_t::load_fmt(RIFF_t &riff)
{
    // WAV file should have a 'fmt ' chunk - throw if fmt is not present
    RIFF_chunk_data_t *fmt = m_fmt(riff.get_root_chunk());
    if (!fmt)
        throw std::runtime_error("File does not have a valid 'fmt ' chunk.");

    memcpy(reinterpret_cast<uint8_t *>(&header), &fmt->get_data().front(), fmt->size());

    // load extra params if present
    if (fmt->get_data().size() > 16)
    {
        // grab size of extra params - magic value 16: end of normal header data
        mempcpy(reinterpret_cast<uint8_t *>(&header.extra_params_size), &fmt->get_data()[16], 2);

        // grab extra params - magic value 18: end of normal header data and size of extra params
        header.extra_params.reserve(header.extra_params_size);
        header.extra_params.insert(header.extra_params.end(), fmt->get_data().begin() + 18, fmt->get_data().end());
    }
}

void WAV_t::load_data(RIFF_t &riff)
{
    // WAV file should have a 'data' chunk - throw if data is not present
    RIFF_chunk_data_t *data = m_data(riff.get_root_chunk());
    if (data == nullptr)
        throw std::runtime_error("File does not have a valid 'data' chunk.");

    std::vector<uint8_t> &d = data->get_data();

    // determine size for sample vector
    // create channels
    samples.insert(samples.begin(), header.num_channels, std::vector<double>());
    // reserve space for samples
    int samples_per_channel{data->size() / header.num_channels};
    for (auto i : samples)
        i.reserve(samples_per_channel);

    // choose encoding type - none is the default
    encoding = WAV_encoding::none;
    if (header.audio_format == FORMAT_PCM) // PCM
    {
        switch (header.bits_per_sample)
        {
        case 8:
            encoding = WAV_encoding::unsigned_8_PCM;
            break;
        case 16:
            encoding = WAV_encoding::signed_16_PCM;
            break;
        case 24:
            encoding = WAV_encoding::signed_24_PCM;
            break;
        case 32:
            encoding = WAV_encoding::signed_32_PCM;
            break;
        }
    }
    else if (header.audio_format == FORMAT_FLOAT) // float
    {
        switch (header.bits_per_sample)
        {
        case 32:
            encoding = WAV_encoding::float_32;
            break;
        case 64:
            encoding = WAV_encoding::float_64;
            break;
        }
    }

    // choose read function based on encoding type
    switch (encoding)
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
        throw std::runtime_error("Unsupported audio encoding format. (" + std::to_string((int)encoding) + ")");
    }
}

void WAV_t::load_cue(RIFF_t &riff)
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

            WAV_t::cue_point point{
                cue_point.sample_start,
                labl_identifiers.find(cue_point.identifier) == labl_identifiers.end() ? std::to_string(cue_point.identifier) : labl_identifiers[cue_point.identifier]};

            cue_points.push_back(point);
        }
    }
}

template <class T>
void WAV_t::load_sample_buffer_int(std::vector<uint8_t> &bytes)
{
    T *buffer = reinterpret_cast<T *>(&bytes.front());
    int channel_counter = 0;
    int total_samples = bytes.size() / sizeof(T);

    // assign each sample to a channel
    for (int i = 0; i < total_samples; i++)
    {
        double new_value = value_map<T, double>(buffer[i], std::numeric_limits<T>::min(), std::numeric_limits<T>::max(), -1.0, 1.0);
        samples[channel_counter++ % header.num_channels].push_back(new_value);
    }
}

void WAV_t::load_sample_buffer_i24(std::vector<uint8_t> &bytes)
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
        double new_value = value_map<int32_t, double>(smp, -0x800000, 0x7fffff, -1.0, 1.0);
        samples[channel_counter++ % header.num_channels].push_back(new_value);
    }
}

template <class T>
void WAV_t::load_sample_buffer_float(std::vector<uint8_t> &bytes)
{
    T *buffer = reinterpret_cast<T *>(&bytes.front());
    int channel_counter = 0;
    int total_samples = bytes.size() / sizeof(T);

    // assign each sample to a channel
    for (int i = 0; i < total_samples; i++)
        samples[channel_counter++ % header.num_channels].push_back(static_cast<double>(buffer[i]));
}

template <class From, class To>
To WAV_t::value_map(From value, From from_min, From from_max, To to_min, To to_max)
{
    double d_value = static_cast<double>(value);
    double d_from_min = static_cast<double>(from_min);
    double d_from_max = static_cast<double>(from_max);
    double d_to_min = static_cast<double>(to_min);
    double d_to_max = static_cast<double>(to_max);

    if (d_value >= d_from_max)
        return to_max;

    if (d_value <= d_from_min)
        return to_min;

    return static_cast<To>((d_value - d_from_min) * (d_to_max - d_to_min) / (d_from_max - d_from_min) + d_to_min);
}

template <class T>
void WAV_t::write_sample_buffer_int(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : samples)
        total_samples += i.size();
    bytes.reserve(total_samples * sizeof(T));

    for (int i = 0; i < total_samples; i++)
    {
        T smp = value_map<double, T>((samples[i % header.num_channels][i / header.num_channels]), -1.0, 1.0, std::numeric_limits<T>::min(), std::numeric_limits<T>::max());

        // push individual bytes from smp into bytes
        for (int j = 0; j < sizeof(T); j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

void WAV_t::write_sample_buffer_i24(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : samples)
        total_samples += i.size();
    bytes.reserve(total_samples * 3);

    for (int i = 0; i < total_samples; i++)
    {
        int32_t smp = value_map<double, int32_t>(samples[i % header.num_channels][i / header.num_channels], -1.0, 1.0, -0x800000, 0x7fffff);

        // push individual bytes from smp into bytes
        for (int j = 0; j < 3; j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

template <class T>
void WAV_t::write_sample_buffer_float(std::vector<uint8_t> &bytes)
{
    // reserve some space
    int total_samples = 0;
    for (auto i : samples)
        total_samples += i.size();
    bytes.reserve(total_samples * sizeof(T));

    for (int i = 0; i < total_samples; i++)
    {
        T smp = static_cast<T>(samples[i % header.num_channels][i / header.num_channels]);

        // push individual bytes from smp into bytes
        for (int j = 0; j < sizeof(T); j++)
            bytes.push_back(reinterpret_cast<uint8_t *>(&smp)[j]);
    }
}

// =========== PUBLIC METHODS ===========

WAV_t::WAV_t() : samples(2, std::vector<double>())
{
    encoding = WAV_encoding::signed_32_PCM;
}

WAV_t::WAV_t(std::string filename)
{
    RIFF_t riff(filename);

    if (strcmp(riff.get_root_chunk().get_form_type(), "WAVE") != 0)
        throw std::runtime_error("File is not a valid WAVE file.");

    load_fmt(riff);
    load_data(riff);
    load_cue(riff);
}

int WAV_t::write(std::string filepath)
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

void WAV_t::update_header()
{
    header.num_channels = samples.size();
    switch (encoding)
    {
    case WAV_encoding::signed_16_PCM:
        header.audio_format = FORMAT_PCM;
        header.bits_per_sample = 16;
        break;
    case WAV_encoding::signed_24_PCM:
        header.audio_format = FORMAT_PCM;
        header.bits_per_sample = 24;
        break;
    case WAV_encoding::signed_32_PCM:
        header.audio_format = FORMAT_PCM;
        header.bits_per_sample = 32;
        break;
    case WAV_encoding::unsigned_8_PCM:
        header.audio_format = FORMAT_PCM;
        header.bits_per_sample = 8;
        break;
    case WAV_encoding::float_32:
        header.audio_format = FORMAT_FLOAT;
        header.bits_per_sample = 32;
        break;
    case WAV_encoding::float_64:
        header.audio_format = FORMAT_FLOAT;
        header.bits_per_sample = 64;
        break;
    }
    calculate_block_align();
    calculate_byte_rate();
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
    switch (encoding)
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
    update_header();
}

std::vector<WAV_t::cue_point> &WAV_t::cues()
{
    return cue_points;
}

// =========== METHODS FOR SAMPLE DATA ===========

uint16_t WAV_t::num_channels() const
{
    return samples.size();
}

std::vector<double> &WAV_t::channel(int i)
{
    return samples[i];
}
std::vector<double> &WAV_t::add_channel()
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
    samples = std::vector<std::vector<double>>(header.num_channels, std::vector<double>());
}

void WAV_t::reset_channel_lengths(double fill)
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