#include <cstdint>
#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include <memory>
#include <exception>

#pragma once

// ====================================================================================================================
/**
 *  A base class for RIFF chunks.
 */
class RIFF_chunk_t
{
protected:
    uint8_t m_identifier[5]{0};
    uint32_t m_size{0};

public:
    virtual ~RIFF_chunk_t() = 0;
    virtual int size() = 0;
    virtual int total_size() = 0;
    virtual void print() = 0;
    virtual void print_full() = 0;
    virtual int write(std::ofstream &f) = 0;

    /**
     * Get the chunk identifier.
     * @return const char * to the chunk identifier.
     */
    const char *get_identifier();

    /**
     * Set the chunk identifier.
     * @param new_id The new identifier for this chunk. An exception will be thrown if an identifier with length != 4 is given.
     */
    void set_identifier(const char *new_id);
};

// ====================================================================================================================
/**
 *  RIFF data chunk class. Chunks containg byte data.
 */
class RIFF_chunk_data_t : public RIFF_chunk_t
{
private:
    std::vector<uint8_t> m_data;

public:
    RIFF_chunk_data_t();
    RIFF_chunk_data_t(const RIFF_chunk_data_t&) = delete;
    RIFF_chunk_data_t(RIFF_chunk_data_t&&) = default;
    
    /**
     * Construct the chunk data from a file.
     * @param f The filestream to read from.
     * @param id Chunk identifier for the new chunk. If nullptr a chunk identifier will be read from the filestream. 
     * An exception will be thrown if an identifier with length != 4 is given.
     */
    RIFF_chunk_data_t(std::ifstream &f, const char *id = nullptr);

    /**
     * Construct a data chunk with given id.
     * @param id The form type for the new chunk. An exception will be thrown if an identifier with length != 4 is given.
     * @param data Optional data to populate the chunk with.
     */
    RIFF_chunk_data_t(const char *id, const std::vector<uint8_t> &data = {});

    ~RIFF_chunk_data_t();

    /**
     * Populate the chunk data from a file.
     * @param f The filestream to read from.
     * @param id Chunk identifier for the new chunk. If nullptr a chunk identifier will be read from the filestream.
     * An exception will be thrown if an identifier with length != 4 is given.
     */
    void read(std::ifstream &f, const char *id);

    /**
     * Write the byte data to the supplied filestream.
     * @param f The filestream to write the bytes to.
     * @return The number of bytes written.
     */
    int write(std::ofstream &f);

    /**
     * Get the currently held chunk data.
     * @return Reference to currently held chunk data.
     */
    std::vector<uint8_t> &get_data();

    /**
     * Set the data for the data chunk.
     * @param new_data The data that replaces the currently held chunk data.
     */
    void set_data(const std::vector<uint8_t> &new_data);

    /**
     * Get the size of the data in the RIFF file in bytes (exluding header information).
     * @see total_size()
     * @return The test results
     */
    int size();

    /**
     * Get the size of the RIFF file in bytes.
     * @see size()
     * @return The test results
     */
    int total_size();

    /**
     * Print basic information about the subchunks.
     */
    void print();

    /**
     * Print information about the subchunks along with complete chunk byte data.
     */
    void print_full();
};

// ====================================================================================================================
/**
 *  RIFF list chunk class. Chunks containing sub chunks (identifier is LIST or RIFF)
 */
class RIFF_chunk_list_t : public RIFF_chunk_t
{
private:
    uint8_t m_form_type[5]{0};
    std::vector<std::unique_ptr<RIFF_chunk_t>> m_subchunks;

public:
    RIFF_chunk_list_t();
    RIFF_chunk_list_t(const RIFF_chunk_list_t&) = delete;
    RIFF_chunk_list_t(RIFF_chunk_list_t&&) = default;

    /**
     * Construct the chunk list from a file.
     * @param f The filestream to read from.
     * @param id Chunk identifier for the new chunk. If nullptr a chunk identifier will be read from the filestream.
     * RIFF specification says list chunks should have "LIST" as the identifier (excluding the root chunk which should have "RIFF" as the identifier).
     * An exception will be thrown if an identifier with length != 4 is given.
     */
    RIFF_chunk_list_t(std::ifstream &f, const char *id = nullptr);

    /**
     * Construct a list chunk with a specific form type.
     * @param form_type The form type for the new chunk. An exception will be thrown if a form type with length != 4 is given.
     */
    RIFF_chunk_list_t(const char *form_type);

    ~RIFF_chunk_list_t();

    /**
     * Get the form type of the chunk.
     * @return const char * to the form type.
     */
    const char *get_form_type();

    /**
     * Set the form type for the chunk.
     * @param new_form_type T
     * @return The test results
     */
    void set_form_type(const char *new_form_type);

    /**
     * Populate the chunk list from a file.
     * @param f The filestream to read from.
     * @param id Chunk identifier for the new chunk. If nullptr a chunk identifier will be read from the filestream.
     * RIFF specification says list chunks should have "LIST" as the identifier (excluding the root chunk which should have "RIFF" as the identifier).
     * An exception will be thrown if an identifier with length != 4 is given.
     */
    void read(std::ifstream &f, const char *id);

    /**
     * Write the byte data to the supplied filestream.
     * @param f The filestream to write the bytes to.
     * @return The number of bytes written.
     */
    int write(std::ofstream &f);

    /**
     * Get a list of the subchunks contained within this LIST chunk.
     * @return A reference to the list of subchunks.
     */
    std::vector<std::unique_ptr<RIFF_chunk_t>> &get_subchunks();

    /**
     * Get the size of the data in the RIFF file in bytes (exluding header information).
     * @see total_size()
     * @return The test results
     */
    int size();

    /**
     * Get the size of the RIFF file in bytes.
     * @see size()
     * @return The test results
     */
    int total_size();

    /**
     * Print basic information about the subchunks.
     */
    void print();

    /**
     * Print information about the subchunks along with complete chunk byte data.
     */
    void print_full();

    /**
     * Tell if a chunk exists within the file structure that matches the specified chunk identifier.
     * @param id The chunk identifier to search for.
     * @return True if the chunk exists, false if it does not.
     */
    bool exists_chunk_with_id(const char *id);

    /**
     * Find the first chunk matching a specified chunk identifier. Pointers returned by this function 
     * may become invalid if significant changes are made to the parent LIST chunk after their creation. 
     * Be careful.
     * @param id The chunk identifier to search for.
     * @return A pointer to the chunk with the specified id. Null pointer if chunk does not exist.
     */
    RIFF_chunk_t *get_chunk_with_id(const char *id);
};

// ====================================================================================================================
/**
 *  RIFF chunk container class. Provides basic operations for accessing chunks within the RIFF file structure.
 */
class RIFF_t
{
private:
    RIFF_chunk_list_t m_riff;
    std::string m_filepath;

public:
    /**
     * Basic constructor. Create an empty RIFF file structure.
     */
    RIFF_t();

    /**
     * Populate a RIFF file structure from a RIFF file.
     * @param filename The RIFF file to parse.
     */
    RIFF_t(std::string filename);

    /**
     * Get the size of the data in the RIFF file in bytes (exluding header information).
     * @see total_size()
     * @return The test results
     */
    int size();

    /**
     * Get the size of the RIFF file in bytes.
     * @see size()
     * @return The test results
     */
    int total_size();

    /**
     * Write the RIFF file to storage. File is written at the currently set file path.
     * @see set_filepath()
     * @return The number of bytes written.
     */
    int write();

    /**
     * Get the current file path of the RIFF file.
     * @return String representation of the location of the RIFF file.
     */
    const std::string &get_filepath();

    /**
     * Set the file path. This has no changes unless write() is called.
     * @param new_file_path New location for the RIFF file.
     * @see write()
     */
    void set_filepath(const std::string &new_file_path);

    /**
     * Print basic information about the RIFF file subchunks.
     */
    void print();

    /**
     * Print information about the RIFF file subchunks along with complete chunk byte data.
     */
    void print_full();

    /**
     * Tell if a chunk exists within the file structure that matches the specified chunk identifier.
     * @param id The chunk identifier to search for.
     * @return True if the chunk exists, false if it does not.
     */
    bool exists_chunk_with_id(const char *id);

    /**
     * Find the first chunk matching a specified chunk identifier. Pointers returned by this function 
     * may become invalid if significant changes are made to the parent LIST chunk after their creation. 
     * Be careful.
     * @param id The chunk identifier to search for.
     * @return A pointer to the chunk with the specified id. Null pointer if chunk does not exist.
     */
    RIFF_chunk_t *get_chunk_with_id(const char *id);

    /**
     * Quick access to the root chunk of the RIFF file.
     * @return Reference to the RIFF root chunk.
     */
    RIFF_chunk_list_t &get_root_chunk();
};