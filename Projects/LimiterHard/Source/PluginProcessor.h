/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "LimiterHard.h"

//==============================================================================
/**
*/
class LimiterHardAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    LimiterHardAudioProcessor();
    ~LimiterHardAudioProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

	int BufferSize;
	float SampleRate;
	int NChannels;
	float LookAhead;
	float PostGain;
	float PreGain;
	void UpdateParameters()
	{
		auto p = Limiter.Algo.GetParameters();
		p.PreGain = PreGain;
		p.PostGain = PostGain;
		Limiter.Algo.SetParameters(p);
	}
	void UpdateCoefficients()
	{
		auto c = Limiter.Algo.GetCoefficients();
		c.LookAheadMS = LookAhead;
		Limiter.Initialize(BufferSize,NChannels, SampleRate, c);
	}
		
private:
	LimiterHardStreaming Limiter;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (LimiterHardAudioProcessor)
};
