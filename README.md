# WAVparser

Read, manipulate, and write wav files in C++.

A general example for reading and writing WAV files using the `WAV_t` class is in the main.cpp file. The sample WAV file included and used in the main.cpp program is from https://freewavesamples.com/ensoniq-zr-76-01-dope-77-bpm.

## Creating a WAV_t object

To create an empty WAV file:

```cpp
WAV_t wav;
```

To parse an existing WAV file:

```cpp
WAV_t wav("path/to/file.wav");
```

To write to a WAV file:

```cpp
wav.write("path/to/new/file.wav");
```

## Header

The header is stored as a private member of the `WAV_t` class.

```cpp
struct WAV_fmt_t
{
    uint16_t audio_format;
    uint16_t num_channels;
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align;
    uint16_t bits_per_sample;
    uint16_t extra_params_size;
    std::vector<uint8_t> extra_params;
};
```

The fields in a `WAV_fmt_t` struct generally depend on things like encoding and sample size. Most of the header details are hidden. Some can be edited:

```cpp
uint16_t sample_rate() const;
void set_sample_rate(uint16_t new_rate);
int sample_size() const;
std::vector<uint8_t> &extra_params(); // extra parameters field in the header, don't know what actually goes here.
```

File encoding can be changed on the fly:

```cpp
WAV_encoding get_encoding() const;
void set_encoding(WAV_encoding new_encoding);
```

block alignment and byte rate can be calculated with some helper functions:

```cpp
int byte_rate = wav.calculate_byte_rate();
int block_align = wav.calculate_block_align();
```

I don't personally see how this information is useful outside of the file header but its available just in case. These functions also set their respective fields in the internal header struct.

## Samples

Individual audio channels and samples can be retreived and assigned using the `channel()` function:

```cpp
int sample = wav.channel(0)[0]; // first sample from channel 0
wav.channel(1)[0] = 0; // first sample from channel 1
```

In stereo audio left channel = 0 and right channel = 1, usually.

To find the number of bytes per sample, use `sample_size()`:

```cpp
int sample_size = wav.sample_size();
```

This will return the number of bytes that the currently set encoding requires per audio sample. Actual audio samples are all stored as `double`s.

To clear all samples from a `WAV_t` object, use `clear_data()`:

```cpp
wav.clear_data();
```

## RIFF_t

WAV files are contained within RIFF files. As such, the `WAV_t` class is built using the `RIFF_t` class. See the [RIFFparser repository]("https://github.com/rami-hansen/RIFFparser") for more information about reading and manipulating RIFF data.