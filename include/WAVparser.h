#include <string>
#include <cstdint>
#include <vector>

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

/**
 * Class for storing and manipulating WAV file data.
 */
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
    /**
     * WAV file header information. Use load_fmt() to load from the RIFF_t
     * object into the header.
     * @see load_fmt();
     */
    WAV_fmt_t header;

    /**
     * Raw samples. Use load_data() to load from the RIFF_t 
     * object into the header.
     * @see load_samples()
     */
    std::vector<uint64_t> samples;

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
     * Load raw byte data from the RIFF_t object into the header.
     */
    void load_fmt();

    /**
     * Load raw byte data from the RIFF_t object into the samples
     * vector.
     */
    void load_data();

    /**
     * Get the raw 'fmt ' data contained in the RIFF_t object. 
     * Changes will not be transferred into the WAV_t object. 
     * Use load_fmt() to transfer from the RIFF_t object. 
     * @return Reference to raw byte data. 
     * @see load_fmt()
     */
    std::vector<uint8_t> &get_fmt();

    /**
     * Get the raw 'data' data contained in the RIFF_t object. 
     * Changes will not be transferred into the WAV_t object. 
     * Use load_data() to transfer from the RIFF_t object. 
     * @return Reference to raw byte data. 
     * @see load_data()
     */
    std::vector<uint8_t> &get_data();

    /**
     * Get the underlying RIFF_t data.
     * @return Reference to the RIFF_t object.
     */
    RIFF_t &get_riff();

    /**
     * @return The number of bytes in a single sample.
     */
    int sample_size();

    /**
     * Get a specific individual sample.
     * @param i The index of the sample to grab
     * @param channel The channel to grab from (default to 0). Requesting a 
     * channel that does not exist will throw an exception. Samples are stored 
     * physically in uint64_t but may be a different logical size. Use 
     * sample_size() to find the true size of a single sample.
     * @see sample_size()
     * @return Reference to the requested sample.
     */
    uint64_t &get_sample(int i, int channel = 0);

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
     * Clear all sample data from the file. Also flushes changes to the RIFF_t object.
     */
    void clear_data();

    /**
     * Set the file path. This has no changes on the underlying data. The 
     * location of the file will be changed when write() is called.
     * @param new_file_path New location for the WAV file.
     * @see write()
     */
    void set_filepath(std::string new_file_path);

    /**
     * Write WAV_t data to the disk at the specified filepath. 
     * Throws exception if no filepath is specified. 
     * @see set_filepath()
     * @return Number of bytes written
     */
    int write();

    /**
     * Quickly print header information
     */
    void print_header();
};