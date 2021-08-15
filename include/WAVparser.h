#include <string>

#include "RIFFparser.h"

#pragma once

class WAV_t
{
private:
    RIFF_t riff;

public:
    WAV_t();
    WAV_t(std::string filename);
};