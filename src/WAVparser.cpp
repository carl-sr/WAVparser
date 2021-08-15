#include "WAVparser.h"

WAV_t::WAV_t()
{
}

WAV_t::WAV_t(std::string filename) : m_riff(filename)
{
}

int WAV_t::write_fmt()
{
}

int WAV_t::write_data()
{
}

void WAV_t::load_fmt()
{
}

void WAV_t::load_data()
{
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