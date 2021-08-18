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
wav.write();
```

Attempting to write a wav file without a specified file path will throw an exception. Set a file path using `set_filepath()`.

## Header

The header is stored as a public member of the `WAV_t` class. Changing header fields without knowing what you're doing can really mess up your audio.

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

// access or assign
int format = wav.header.audio_format;
wav.header.num_channels = 1;
```

block alignment and byte rate can be calculated with some helper functions:

```cpp
int byte_rate = wav.calculate_byte_rate();
int block_align = wav.calculate_block_align();
```

These helper functions also set their respective fields in the header.

Raw header 'fmt ' data can be accessed using the `get_fmt()` function. Changes made here will not be reflected in the `WAV_t` header unless `load_fmt()` is called.

```cpp
std::vector<uint8_t> &raw_header_data = get_fmt();
raw_header_data[0] = 0xff;
wav.load_fmt();
```

## Samples

Individual audio samples can be retreived and assigned using the `get_sample()` function:

```cpp
int sample = wav.get_sample(0, 0); // first sample from channel 0
wav.get_sample(0, 1) = 0; // first sample from channel 1
```

The size of a sample is stored internally as a `uint64_t`. Actual size varies depending on the wav file. Sample sizes do not exceed 64 bits, as far as I know...

To find the number of bytes per sample, use `sample_size()`:

```cpp
int sample_size = wav.sample_size();
```

Access to raw 'data' section bytes is gained with `get_data()`. Like getting raw header data, changes made here will not affect the samples vector unless `load_data()` is called.

```cpp
std::vector<uint8_t> &raw_data = get_data();
raw_data[0] = 0xff;
wav.load_data();
```

To clear all samples from a `WAV_t` object, use `clear_data()`:

```cpp
wav.clear_data();
```

`clear_data()` clears both the samples vector as well as the underlying byte data held in the RIFF_t object.

## RIFF_t

WAV files are contained within RIFF files. As such, the `WAV_t` class is built using the `RIFF_t` class. Access the underlying `RIFF_t` object using the `get_riff()` method. See the
[RIFFparser repository](https://github.com/rami-hansen/RIFFparser)
for more information about reading and manipulating RIFF data.
