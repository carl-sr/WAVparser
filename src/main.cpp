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

    // get the total number of samples (from all channels)
    int num_samples{0};
    for (int i = 0; i < wav.num_channels(); i++)
        num_samples += wav.channel(i).size();

    // make sure each channel has an equal number of samples
    wav.reset_channel_lengths();
    
    // remove half of the samples
    int samples_per_channel = wav.channel(0).size();
    for (int i = 0; i < wav.num_channels(); i++)
    {
        std::vector<float> &channel = wav.channel(i);
        channel.erase(channel.begin(), channel.begin() + samples_per_channel / 2);
    }

    // get a single sample (first sample in channel 0)
    float sample = wav.channel(0)[0];

    // set a single sample (first sample in channel 0)
    wav.channel(0)[0] = 0.0f;

    // write the WAV file to disk
    wav.write("half.wav");

    // clear the wav file sample data
    wav.clear_data();

    wav.write("empty.wav");

    return 0;
}