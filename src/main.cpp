#include <iostream>

#include "./WAVparser.h"

int main(int argc, char *argv[])
{
    std::cout << "WAVparser" << std::endl;
    
    // Create a WAV_t object from the file 'sample.wav'
    // sample.wav from https://freewavesamples.com/ensoniq-zr-76-01-dope-77-bpm
    WAV_t wav("./sample.wav");

    // quick print header information
    wav.print_header();

    // vector of individual samples is available as the public member samples
    int sample_count = wav.samples.size();

    // header information is also available as a public member of type WAV_fmt_t
    int byte_rate = wav.header.sample_rate;

    // remove half of the samples
    wav.samples.erase(wav.samples.begin() + (sample_count / 2), wav.samples.end());

    // get a single sample (first sample in channel 0)
    uint64_t smp = wav.get_sample(0, 0);

    // set a single sample (first sample in channel 0)
    wav.get_sample(0) = 0;

    // set a new filepath
    wav.set_filepath("./new_wav.wav");

    // write the WAV file to disk
    wav.write();

    // access is also provided to the underlying RIFF structure
    // https://github.com/rami-hansen/RIFFparser
    RIFF_t &riff = wav.get_riff();

    // print riff information
    riff.print();

    // get raw fmt and data
    std::vector<uint8_t> &fmt = wav.get_fmt();
    std::vector<uint8_t> &data = wav.get_data();

    // clear the wav file sample data
    wav.clear_data();

    return 0;
}