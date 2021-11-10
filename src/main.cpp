#include <iostream>

#include "./WAVparser.h"

int main(int argc, char *argv[])
{
    std::cout << "WAVparser" << std::endl;

    // Create a WAV_t object from the file 'sample.wav'
    // sample.wav from https://freewavesamples.com/ensoniq-zr-76-01-dope-77-bpm
    WAV_t wav("./sample.wav");

    // print all of the samples in the file
    for(int i = 0; i < wav.channel(0).size(); i++)
    {
        for(int j = 0; j < wav.get_num_channels(); j++)
        {
            double smp = wav.channel(j)[i];
            printf(smp >= 0.0f ? " %.8f  " : "%.8f  ", smp);
        }
        putchar('\n');
    }
    
    // quick print header information
    wav.print_header();

    // set encodings and write files
    wav.set_encoding(WAV_encoding::unsigned_8_PCM);
    wav.write("out_u8.wav");

    wav.set_encoding(WAV_encoding::signed_16_PCM);
    wav.write("out_i16.wav");

    wav.set_encoding(WAV_encoding::signed_24_PCM);
    wav.write("out_i24.wav");

    wav.set_encoding(WAV_encoding::signed_32_PCM);
    wav.write("out_i32.wav");

    wav.set_encoding(WAV_encoding::float_32);
    wav.write("out_f32.wav");

    wav.set_encoding(WAV_encoding::float_64);
    wav.write("out_f64.wav");


    // get the total number of samples (from all channels)
    int num_samples{0};
    for (int i = 0; i < wav.get_num_channels(); i++)
        num_samples += wav.channel(i).size();

    // make sure each channel has an equal number of samples
    wav.reset_channel_lengths();

    // remove half of the samples
    int samples_per_channel = wav.channel(0).size();
    for (int i = 0; i < wav.get_num_channels(); i++)
    {
        std::vector<double> &channel = wav.channel(i);
        channel.erase(channel.begin(), channel.begin() + samples_per_channel / 2);
    }

    // get a single sample (first sample in channel 0)
    double sample = wav.channel(0)[0];

    // set a single sample (first sample in channel 0)
    wav.channel(0)[0] = 0.0f;

    // write the WAV file to disk
    wav.write("sample_half.wav");

    // clear the wav file sample data
    wav.clear_data();

    wav.write("sample_empty.wav");

    return 0;
}