#include <iostream>

#include "./WAVparser.h"

int main()
{
    std::cout << "WAVparser" << std::endl;
    // ================================================================================================
    // ================== CREATING WAV OBJECTS ==================
    // Create a WAV_t object from the file 'sample.wav'
    // sample.wav from https://freewavesamples.com/ensoniq-zr-76-01-dope-77-bpm
    // different sample types are supported - stick to float and double, other's don't work very well.
    WAV<float> floatWav("./sample.wav");
    WAV<double> doubleWav("./sample.wav");

    // ================================================================================================
    // ================== HEADER METHODS ==================
    // quick print header information
    printf("Header information:\n");
    floatWav.printHeader();

    // get access to the raw header data. Generally a bad idea.
    // Messing with this without knowing what you're doing can
    // result in a corrupted WAV file.
    WAV_fmt_t floatWavHeader = floatWav.getRawHeader();

    // To get the WAV object's sample rate:
    floatWav.getSampleRate();

    // To set the WAV ojects's sample rate:
    floatWav.setSampleRate(44100);

    // setSampleRate only sets the header value. To stretch the WAV
    // so that it sounds similar but with a different sample rate
    // use stretchToSampleRate() (INCOMPLETE):
    // floatWav.stretchToSampleRate(44100);

    // To get the number of bytes an encoded sample will take:
    floatWav.getSampleSize();

    // Some WAV files hold extra information in the header in a field
    // called extra params. To access this information:
    floatWav.getExtraParams();

    // ================================================================================================
    // ================== ENCODINGS ==================

    // The supported encoding types are:
    //   - Unsigned 8 bit PCM
    //   - Signed 16 bit PCM
    //   - Signed 24 bit PCM
    //   - Signed 32 bit PCM
    //   - 32 bit floating point
    //   - 64 bit floating point

    // To get the current objects encoding:
    WAV_encoding encoding = floatWav.getEncoding();

    // To get a list of supported encoding types and their
    // respective enum values, use getAvailableEncodings():
    printf("\nEncodings:\n");
    for (auto encoding : floatWav.getAvailableEncodings())
        printf("%s = %i\n", encoding.encoding_str, (int)encoding.encoding_type);

    // Encodings can be set using setEncoding()
    // This code writes a file using each of the available encoding:
    for (auto encoding : floatWav.getAvailableEncodings())
    {
        floatWav.setEncoding(encoding.encoding_type);
        floatWav.write(std::string(encoding.encoding_str) + ".wav");
    }

    // ================================================================================================
    // ================== CUE POINTS ==================
    // If your WAV file has cue points they are accessible and editable
    // using getCuePoints():
    auto cues = floatWav.getCuePoints();

    // Any cue points added to the object will be included in the file
    // written to disk.
    WAV<float>::cue_point newCuePoint;
    newCuePoint.sample_offset = 0;
    newCuePoint.label = "Beginning of file";
    cues.push_back(newCuePoint);

    // ================================================================================================
    // ================== MANIPULATING CHANNELS ==================
    // The number of channels can be increased and decreased.
    // Get the number of channels:
    int numChannels = floatWav.getNumChannels();

    printf("\nBefore channel operations numChannels: %i\n", doubleWav.getNumChannels());
    // Remove a channel:
    doubleWav.removeChannel(0);
    printf("After removing a channel numChannels: %i\n", doubleWav.getNumChannels());

    // Set the number of channels:
    floatWav.setNumChannels(2);
    printf("After setting numChannels to 2 numChannels: %i\n", doubleWav.getNumChannels());

    // Added channels are initialized to zero
    // Add a channel:
    doubleWav.addChannel();
    printf("After adding a channel numChannels: %i\n", doubleWav.getNumChannels());


    // The object can be converted to mono. This averages samples in each channel
    // to create a single channel that still sounds the same:
    // floatWav.convertToMono();

    // To empty all samples from the object:
    doubleWav.clearData();

    // ================================================================================================
    // ================== MANIPULATING SAMPLES ==================

    // get a specific channel using the sampleAt function
    static const int channel = 0;
    static const int sampleIndex = 100;
    float sample = floatWav.sampleAt(channel, sampleIndex);

    // A good way to manipulate WAV objects is by using iterators
    // Two types of iterators are supplied:
    //   - Cross channel iterators
    //   - Single channel iterators

    // Cross channel iterators
    WAV<float>::iterator iterBegin = floatWav.begin();
    WAV<float>::iterator iterEnd = floatWav.end();

    // Individual channels at each sample index can be accessed
    // using the channel() method:
    float leftChannelFirstSample = iterBegin->channel(0);
    float rightChennelFirstSample = iterBegin->channel(1);

    // Some of the usual iterator operations apply:
    WAV<float> newWav;
    newWav.insert(newWav.begin(), floatWav.begin(), floatWav.begin() + 1000);
    newWav.erase(newWav.begin(), newWav.begin() + 980);
    newWav.write("iterator.wav"); // this will have 20 audio samples

    // print all samples using iterators:
    printf("\nSamples inserted using iterators:\n");
    for (auto iter : newWav)
    {
        for (int ch = 0; ch < floatWav.getNumChannels(); ch++)
            printf("%s%f ", iter.channel(ch) < 0.0f ? "" : " ", iter.channel(ch));
        printf("\n");
    }

    // Single channel iterators
    WAV<float>::channelIterator leftIter = newWav.channelBegin(0);

    printf("\nRight channel samples:\n");
    for(auto rightIter = newWav.channelBegin(1); rightIter != newWav.channelEnd(1); rightIter++)
        printf("%f\n", *rightIter);

    printf("\n");

    return 0;
}