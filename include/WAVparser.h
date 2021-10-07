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
class WAV_t
{
public:
    struct cue_point
    {
        uint32_t sample_offset{0};
        std::string label;
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
    To value_map(From value, From from_min, From from_max, To to_min, To to_max);

    std::vector<std::vector<double>> samples;
    WAV_fmt_t header;
    WAV_encoding encoding{WAV_encoding::none};

    std::vector<cue_point> cue_points;

public:
    /**
     * Construct an empty WAV file. Contains no samples and 
     * default header data.
     */
    WAV_t();

    /**
     * Construct a WAV_t object from a WAV file.
     * @param filename The file to parse.
     */
    WAV_t(std::string filename);

    /**
     * Write WAV_t data to the disk at the specified filepath. 
     * @return Number of bytes written
     */
    int write(std::string filepath);

    // =========== METHODS FOR HEADER INFORMATION ===========

    /**
     * Get the current sample rate of the WAV data.
     */
    uint16_t get_sample_rate() const;

    /**
     * Set a new sample rate for the WAV data.
     * @param new_rate The new sample rate.
     */
    void set_sample_rate(uint16_t new_rate);

    /**
     * Get the size in bytes of a sample when written.
     */
    int sample_size() const;

    /**
     * Get the vector containing byte data for extra parameters. 
     * Extra parameters are not parsed by the object.
     */
    std::vector<uint8_t> &extra_params();

    /**
     * Quickly print header information
     */
    void print_header();

    /**
     * Helper function to set header byte rate. Generally only used internally but can be useful to do 
     * correct calculations when changing header information manually.
     * @return The new calculated byte rate.
     */
    uint32_t calculate_byte_rate();

    /**
     * Helper function to set header block align. Generally only used internally but can be useful to do 
     * correct calculations when changing header information manually.
     * @return The new calculated block align.
     */
    uint16_t calculate_block_align();

    /**
     * Get the current encoding type of the WAV data.
     */
    WAV_encoding get_encoding() const;

    /**
     * Set a new encoding type for the WAV data
     * @param new_encoding The encoding type to use during the next write.
     */
    void set_encoding(WAV_encoding new_encoding);

    /**
     * Get the list of cue points.
     * @return Reference to list of cue points.
     */
    std::vector<WAV_t::cue_point> &cues();

    // =========== METHODS FOR SAMPLE DATA ===========

    /**
     * Get the current number of audio channels.
     */
    uint16_t get_num_channels() const;

    /**
     * Set the number of channels
     * If the new number of channels is smaller than the current number, channels will be deleted.
     * If it is greater than the current number, new channels initialized to zero will be created.
     */
    void set_num_channels(int num_channels);

    /**
     * Get the current number of samples in the sample buffer
     */
    int get_num_samples();

    /**
     * Get channel data.
     * @param i The channel number to grab
     */
    std::vector<double> &channel(int i);

    /**
     * Add an audio channel.
     * @returns Reference to the newly created audio channel
     */
    std::vector<double> &add_channel();

    /**
     * Remove an audio channel.
     * @param i The index of the channel to remove
     */
    void remove_channel(int i);

    /**
     * Swap two audio channels.
     * @param a First channel to swap.
     * @param b Second channel to swap.
     */
    void swap_channels(int a, int b);

    /**
     * Clear all sample data from the file.
     */
    void clear_data();

    /**
     * Set the length of all channels to the length of the longest channel.
     * @param fill The value to append to short channels.
     */
    void reset_channel_lengths(double fill = 0.0f);
};