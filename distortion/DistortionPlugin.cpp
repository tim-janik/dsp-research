#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_audio_utils/juce_audio_utils.h>

#include "distortiondsp.hh"

using namespace juce; // This allows you to use 'AudioProcessor' instead of 'juce::AudioProcessor'
//==============================================================================
class DistortionProcessor final : public AudioProcessor
{
  DistortionDSP distortion_dsp;
public:
    //==============================================================================
    DistortionProcessor()
        : AudioProcessor (BusesProperties().withInput  ("Input",  AudioChannelSet::stereo())
                                           .withOutput ("Output", AudioChannelSet::stereo()))
    {
      addParameter (drive = new AudioParameterFloat ({ "drive", 1 }, "Drive", -6.0f, 36.0f, 0.0f));
      addParameter (symmetry = new AudioParameterFloat ({ "symmetry", 1 }, "Symmetry", -100.0f, 100.0f, 0.0f));
      addParameter (mode = new AudioParameterChoice ({ "mode", 1 }, "Distortion Mode",
        {
          "tanh",
          "sin",
          "west-coast+adaa",
          "west-coast+adaa+lpf",
          "hard-clip",
          "soft-clip3",
          "soft-clip4",
          "soft-clip5"
        }, 0));
      addParameter (oversample_param = new AudioParameterChoice ({ "oversample", 1 }, "Oversample", { "1x", "2x", "4x", "8x" }, 0));
      addParameter (slew = new AudioParameterFloat ({ "slew", 1 }, "Slew", 0.0f, 100.0f, 100.0f));

      auto freq_range = NormalisableRange<float>(
        20.0f,
        20000.0f,
        [](float start, float end, float proportion)
          {
            return start * std::pow (end / start, proportion);
          },
        [](float start, float end, float value)
          {
            return std::log (value / start) / std::log (end / start);
          }
      );
      auto lp_freq_range = NormalisableRange<float>(
        500.0f,
        20000.0f,
        [](float start, float end, float proportion)
          {
            return start * std::pow (end / start, proportion);
          },
        [](float start, float end, float value)
          {
            return std::log (value / start) / std::log (end / start);
          }
      );
      auto hp_freq_range = NormalisableRange<float>(
        20.0f,
        2000.0f,
        [](float start, float end, float proportion)
          {
            return start * std::pow (end / start, proportion);
          },
        [](float start, float end, float value)
          {
            return std::log (value / start) / std::log (end / start);
          }
      );
      auto Q_range = NormalisableRange<float>(
        0.1f,
        10.0f,
        [](float start, float end, float proportion)
          {
            return start * std::pow (end / start, proportion);
          },
        [](float start, float end, float value)
          {
            return std::log (value / start) / std::log (end / start);
          }
      );

      addParameter (pre_eq_freq = new AudioParameterFloat ({ "pre_eq_freq", 1 }, "EQ Freq", freq_range, 1000.0f));
      addParameter (pre_eq_gain = new AudioParameterFloat ({ "pre_eq_gain", 1 }, "EQ Gain", -24.f, 24.f, 12.f));
      addParameter (pre_eq_Q = new AudioParameterFloat ({"pre_eq_Q", 1 }, "EQ Q", Q_range, 1.f));
      addParameter (post_hp_freq = new AudioParameterFloat ({"post_hp_freq", 1 }, "Post HP Freq", hp_freq_range, 80.f));
      addParameter (post_lp_freq = new AudioParameterFloat ({"post_lp_freq", 1 }, "Post LP Freq", lp_freq_range, 8000.f));

      addParameter (width = new AudioParameterFloat ({ "width", 1 }, "Width", 0.f, 150.f, 100.0f));
      addParameter (wet_gain = new AudioParameterFloat ({ "wet_gain", 1 }, "Wet Gain", -24.0f, 24.0f, 0.0f));
      addParameter (mix = new AudioParameterFloat ({ "mix", 1 }, "Mix", 0.0f, 100.0f, 100.0f));
    }

    //==============================================================================
    void prepareToPlay (double sample_rate, int) override
    {
      distortion_dsp.reset (sample_rate);
      updateParams (true);
    }
    void releaseResources() override {}

   void
   updateParams (bool now)
   {
     int new_oversample = 1;
     if (oversample_param->getIndex() == 0)
       new_oversample = 1;
     if (oversample_param->getIndex() == 1)
       new_oversample = 2;
     if (oversample_param->getIndex() == 2)
       new_oversample = 4;
     if (oversample_param->getIndex() == 3)
       new_oversample = 8;
     distortion_dsp.set_oversample (new_oversample);
     distortion_dsp.set_drive (drive->get(), now);
     distortion_dsp.set_symmetry (symmetry->get(), now);
     distortion_dsp.set_slew (slew->get(), now);
     distortion_dsp.set_mode (mode->getIndex());
     distortion_dsp.set_pre_eq_params (pre_eq_freq->get(), pre_eq_gain->get(), pre_eq_Q->get(), now);
     distortion_dsp.set_post_lp (post_lp_freq->get(), now);
     distortion_dsp.set_post_hp (post_hp_freq->get(), now);
     distortion_dsp.set_width (width->get(), now);
     distortion_dsp.set_wet_gain (wet_gain->get(), now);
     distortion_dsp.set_mix (mix->get(), now);
   }

    void processBlock (AudioBuffer<float>& buffer, MidiBuffer&) override
    {
      updateParams (false);

      float *left = buffer.getWritePointer (0);
      float *right = buffer.getWritePointer (1);

      distortion_dsp.process_block (left, right, buffer.getNumSamples());
    }

    bool supportsDoublePrecisionProcessing() const override
    {
      return false;
    }
    void processBlock (AudioBuffer<double>& buffer, MidiBuffer&) override
    {
    }

    //==============================================================================
    AudioProcessorEditor* createEditor() override          { return new GenericAudioProcessorEditor (*this); }
    bool hasEditor() const override                        { return true;   }

    //==============================================================================
    const String getName() const override                  { return "Distortion Plugin"; }
    bool acceptsMidi() const override                      { return false; }
    bool producesMidi() const override                     { return false; }
    double getTailLengthSeconds() const override           { return 0; }

    //==============================================================================
    int getNumPrograms() override                          { return 1; }
    int getCurrentProgram() override                       { return 0; }
    void setCurrentProgram (int) override                  {}
    const String getProgramName (int) override             { return "None"; }
    void changeProgramName (int, const String&) override   {}

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override
    {
        //MemoryOutputStream (destData, true).writeFloat (*gain);
    }

    void setStateInformation (const void* data, int sizeInBytes) override
    {
        //gain->setValueNotifyingHost (MemoryInputStream (data, static_cast<size_t> (sizeInBytes), false).readFloat());
    }

    //==============================================================================
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override
    {
        const auto& mainInLayout  = layouts.getChannelSet (true,  0);
        const auto& mainOutLayout = layouts.getChannelSet (false, 0);

        return (mainInLayout == mainOutLayout && (! mainInLayout.isDisabled()));
    }

private:
    //==============================================================================
    AudioParameterFloat* drive;
    AudioParameterFloat* symmetry;
    AudioParameterFloat* mix;
    AudioParameterFloat* width;
    AudioParameterFloat* wet_gain;
    AudioParameterFloat* slew;
    AudioParameterFloat* pre_eq_freq;
    AudioParameterFloat* pre_eq_gain;
    AudioParameterFloat* pre_eq_Q;
    AudioParameterFloat* post_lp_freq;
    AudioParameterFloat* post_hp_freq;
    AudioParameterChoice* mode;
    AudioParameterChoice* oversample_param;

#if 0
    DistortionDSP saturation;
#endif

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (DistortionProcessor)
};

// This must be outside the class definition
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new DistortionProcessor();
}
