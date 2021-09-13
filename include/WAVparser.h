#include <string>
#include <cstdint>
#include <vector>
#include <limits>
#include <iostream>

#include "RIFFparser.h"

#pragma once

#pragma pack(2)

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
private:
    // quick access
    RIFF_chunk_data_t *m_data(RIFF_chunk_list_t &data);
    RIFF_chunk_data_t *m_fmt(RIFF_chunk_list_t &fmt);

    // write fmt and data sections into m_riff
    int write_fmt(RIFF_t &riff);
    int write_data(RIFF_t &riff);

    // Load raw byte data from the RIFF_t object into the header
    void load_fmt(RIFF_chunk_data_t &riff);

    // Load raw byte data from the RIFF_t object into the samples vector.
    void load_data(RIFF_chunk_data_t &riff);

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
    uint16_t sample_rate() const;

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
     * Set the size in bytes of a sample when written.
     * @param new_size The new size of a sample.
     */
    void set_sample_size(int new_size);

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

    // =========== METHODS FOR SAMPLE DATA ===========

    /**
     * Get the current number of audio channels.
     */
    uint16_t num_channels() const;

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