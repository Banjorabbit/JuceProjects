/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PitchShiftAdaptiveResolution.h"

//==============================================================================
/**
*/
class PitchShiftAdaptiveResolutionPluginAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    PitchShiftAdaptiveResolutionPluginAudioProcessor();
    ~PitchShiftAdaptiveResolutionPluginAudioProcessor();

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

	float Pitch, DetectTransient;
	void UpdateParameters()
	{
		auto p = PitchShifter0.Algo.GetParameters();
		p.StretchFactor = Pitch;
		PitchShifter0.Algo.SetParameters(p);
		PitchShifter1.Algo.SetParameters(p);
		auto p2 = PitchShifter0.Algo.TransientDetection.GetParameters();
		p2.TransientThreshold = DetectTransient;
		PitchShifter0.Algo.TransientDetection.SetParameters(p2);
		PitchShifter1.Algo.TransientDetection.SetParameters(p2);
	}
private:
	PitchShiftAdaptiveResolutionStreaming PitchShifter0;
	PitchShiftAdaptiveResolutionStreaming PitchShifter1;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PitchShiftAdaptiveResolutionPluginAudioProcessor)
};
