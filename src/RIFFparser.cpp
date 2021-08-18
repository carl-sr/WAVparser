#include "RIFFparser.h"

// ====================================================================================================================
RIFF_chunk_t::~RIFF_chunk_t() {}

void RIFF_chunk_t::set_identifier(const char *new_id)
{
    if (strlen(new_id) != 4)
        throw std::invalid_argument("RIFF chunk identifier must be exactly four characters.");

    m_identifier[0] = new_id[0];
    m_identifier[1] = new_id[1];
    m_identifier[2] = new_id[2];
    m_identifier[3] = new_id[3];
}

const char *RIFF_chunk_t::get_identifier()
{
    return reinterpret_cast<const char *>(m_identifier);
}

// ====================================================================================================================
RIFF_chunk_data_t::RIFF_chunk_data_t(std::ifstream &f, const char *id)
{
    read(f, id);
}

RIFF_chunk_data_t::RIFF_chunk_data_t()
{
}

RIFF_chunk_data_t::RIFF_chunk_data_t(const char *id, const std::vector<uint8_t> &data)
{
    set_identifier(id);
    set_data(data);
}

RIFF_chunk_data_t::~RIFF_chunk_data_t()
{
}

void RIFF_chunk_data_t::read(std::ifstream &f, const char *id)
{
    // if no id is supplied, read it from the filestream
    if (id)
    {
        set_identifier(id);
    }
    else
    {
        char identifier[5]{0};
        f.read(identifier, 4);
        set_identifier(identifier);

    }
    f.read(reinterpret_cast<char *>(&m_size), sizeof(m_size));

    // grab data
    m_data.reserve(m_size);
    for (int i = 0; i < m_size; i++)
        m_data.push_back(f.get());

    if(m_data.size() > 64)
    {
        putchar('\n');
        return;
    }
}

int RIFF_chunk_data_t::write(std::ofstream &f)
{
    f.write(reinterpret_cast<const char *>(m_identifier), 4);
    m_size = m_data.size();
    f.write(reinterpret_cast<const char *>(&m_size), 4);
    for (auto i : m_data)
        f.write(reinterpret_cast<const char *>(&i), 1);

    int bytes = 8 + m_data.size();

    // padding byte if data is odd sized
    if (bytes % 2 != 0)
    {
        const char pad{'\0'};
        f.write(&pad, 1);
        bytes++;
    }
    return bytes;
}

void RIFF_chunk_data_t::set_data(const std::vector<uint8_t> &new_data)
{
    m_data = new_data;
}

std::vector<uint8_t> &RIFF_chunk_data_t::get_data()
{
    return m_data;
}

int RIFF_chunk_data_t::size()
{
    return m_data.size();
}

int RIFF_chunk_data_t::total_size()
{
    int bytes = m_data.size() + 8;

    // padding byte
    if (bytes % 2 != 0)
        bytes += 1;

    return bytes;
}

void RIFF_chunk_data_t::print()
{
    printf("RIFF_chunk_data: (length %lu) id: %s\n", m_data.size(), get_identifier());
    for (int i = 0; i < m_data.size() && i < 8; i++)
    {
        printf(" %02x", m_data[i]);
    }

    if (m_data.size() > 8)
        printf(" ... ");

    putchar('\n');
}

void RIFF_chunk_data_t::print_full()
{
    printf("RIFF_chunk_data: (length %lu) id: %s\n", m_data.size(), get_identifier());
    for (int i = 0; i < m_data.size(); i++)
    {
        if (i % 8 == 0 && i != 0)
            putchar('\n');

        printf(" %02x", m_data[i]);
    }
    putchar('\n');
}

// ====================================================================================================================
RIFF_chunk_list_t::RIFF_chunk_list_t()
{
}

RIFF_chunk_list_t::RIFF_chunk_list_t(std::ifstream &f, const char *id)
{
    read(f, id);
}

RIFF_chunk_list_t::RIFF_chunk_list_t(const char *form_type)
{
    set_identifier("LIST");
    set_form_type(form_type);
}

RIFF_chunk_list_t::~RIFF_chunk_list_t()
{
}

const char *RIFF_chunk_list_t::get_form_type()
{
    return reinterpret_cast<const char *>(m_form_type);
}

void RIFF_chunk_list_t::set_form_type(const char *new_form_type)
{
    if (strlen(new_form_type) != 4)
        throw std::runtime_error("RIFF chunk form type must be exactly four characters.");

    m_form_type[0] = new_form_type[0];
    m_form_type[1] = new_form_type[1];
    m_form_type[2] = new_form_type[2];
    m_form_type[3] = new_form_type[3];
}

void RIFF_chunk_list_t::read(std::ifstream &f, const char *id)
{
    if (id)
        set_identifier(id);

    else
    {
        char identifier[5]{0};
        f.read(identifier, 4);
        set_identifier(identifier);
    }

    f.read(reinterpret_cast<char *>(&m_size), sizeof(m_size));
    f.read(reinterpret_cast<char *>(&m_form_type), 4);

    // size == 4 means the list contains only the form type
    // form type has already been read, return
    if(m_size <= 4)
        return;

    // while bytes from original position until position + chunk size have been read, or file has been read completely
    std::streampos fpos = f.tellg();
    while (!f.eof() && f.tellg() <= fpos + static_cast<std::streampos>(m_size))
    {
        char identifier[5]{0};

        // skip padding bytes, eof sometimes doesn't work in the initial while test? check it here also.
        if(f.peek() == 0 || f.eof())
        {
            f.read(identifier, 1);
            continue;
        }

        f.read(identifier, 4);

        // determine which type of chunk to add based on its identifier
        if (strcmp(identifier, "LIST") == 0)
            m_subchunks.push_back(std::make_unique<RIFF_chunk_list_t>(f, identifier));
        else
            m_subchunks.push_back(std::make_unique<RIFF_chunk_data_t>(f, identifier));
    }
}

int RIFF_chunk_list_t::write(std::ofstream &f)
{
    int bytes{0};
    f.write(reinterpret_cast<const char *>(m_identifier), 4);

    // size of contained data minus size of header bytes
    uint32_t size = total_size() - 8;

    // in case a padding byte is needed at the end
    if (size % 2 != 0)
        size += 1;

    f.write(reinterpret_cast<const char *>(&size), 4);
    f.write(reinterpret_cast<const char *>(&m_form_type), 4);

    bytes += 12;

    for (auto &i : m_subchunks)
        bytes += i.get()->write(f);

    // if the calculated size of the bytes does not match - write pad bytes until it does
    for (int i = 0; i + bytes < size + 8; i++)
    {
        const char pad{'\0'};
        f.write(&pad, 1);
        bytes++;
    }

    return bytes;
}

std::vector<std::unique_ptr<RIFF_chunk_t>> &RIFF_chunk_list_t::get_subchunks()
{
    return m_subchunks;
}

int RIFF_chunk_list_t::size()
{
    int bytes{0};
    for (auto &i : m_subchunks)
        bytes += i.get()->size();

    return bytes;
}

int RIFF_chunk_list_t::total_size()
{
    int bytes{0};
    for (auto &i : m_subchunks)
        bytes += i.get()->total_size();

    // compensate for padding bytes if the number of bytes is odd
    if (bytes % 2 != 0)
        bytes += 1;

    // include size of id, size, form type
    return bytes + 12;
}

void RIFF_chunk_list_t::print()
{
    printf("\nRIFF_chunk_list: (length %lu chunks) id: %s, form type: %s\n", m_subchunks.size(), get_identifier(), get_form_type());
    for (auto &i : m_subchunks)
    {
        i.get()->print();
    }
}

void RIFF_chunk_list_t::print_full()
{
    printf("\nRIFF_chunk_list: (length %lu chunks) id: %s, form type: %s\n", m_subchunks.size(), get_identifier(), get_form_type());
    for (auto &i : m_subchunks)
    {
        i.get()->print_full();
    }
}

bool RIFF_chunk_list_t::exists_chunk_with_id(const char *id)
{
    return get_chunk_with_id(id) ? true : false;
}

RIFF_chunk_t *RIFF_chunk_list_t::get_chunk_with_id(const char *id)
{
    if (strlen(id) != 4)
        return nullptr;

    for (auto &i : m_subchunks)
    {
        // compare identifiers
        RIFF_chunk_t *current = i.get();
        if (strcmp(current->get_identifier(), id) == 0)
            return current;

        // if the current chunk is a list chunk, recurse into it
        RIFF_chunk_list_t *dyn_next = dynamic_cast<RIFF_chunk_list_t *>(current);
        if (dyn_next)
        {
            RIFF_chunk_t *next = dyn_next->get_chunk_with_id(id);
            if (next)
                return next;
        }
    }
    return nullptr;
}

// ====================================================================================================================
RIFF_t::RIFF_t()
{
    m_riff.set_identifier("RIFF");
    m_riff.set_form_type("NULL");
}

RIFF_t::RIFF_t(std::string filename)
{
    m_filepath = filename;
    
    std::ifstream f(filename, std::ios::binary);
    if (!f.is_open())
        throw std::runtime_error("An error occurred opening the specified RIFF file.");

    // read file identifier, verify that this is a valid RIFF file
    char identifier[5]{0};
    f.read(identifier, 4);

    if (strcmp(identifier, "RIFF") != 0)
        throw std::runtime_error("The specified file is not a valid RIFF file.");

    // this is a valid riff file, read the rest of it
    m_riff.read(f, identifier);

    f.close();
}

int RIFF_t::size()
{
    return m_riff.size();
}

int RIFF_t::total_size()
{
    int size = m_riff.total_size();

    // compensate for padding bytes if the number of bytes is odd
    if (size % 2)
        size += 1;

    return size;
}

void RIFF_t::print()
{
    m_riff.print();
}

void RIFF_t::print_full()
{
    m_riff.print_full();
}

int RIFF_t::write()
{
    if (m_filepath.empty())
        throw std::runtime_error("No file path specified.");

    int bytes{0};

    std::ofstream f(m_filepath, std::ios::binary | std::ios::trunc);

    if (!f.is_open())
        throw std::runtime_error("Unable to open specified file for writing.");

    bytes += m_riff.write(f);

    return bytes;
}

const std::string &RIFF_t::get_filepath()
{
    return m_filepath;
}

void RIFF_t::set_filepath(const std::string &new_file_path)
{
    m_filepath = new_file_path;
}

bool RIFF_t::exists_chunk_with_id(const char *id)
{
    if (strlen(id) != 4)
        return false;

    return get_chunk_with_id(id) ? true : false;
}

RIFF_chunk_t *RIFF_t::get_chunk_with_id(const char *id)
{
    if (strlen(id) != 4)
        return nullptr;

    // if the root RIFF chunk is requested
    if (strcmp(id, "RIFF") == 0)
        return &m_riff;

    return m_riff.get_chunk_with_id(id);
}

RIFF_chunk_list_t &RIFF_t::get_root_chunk()
{
    return m_riff;
}