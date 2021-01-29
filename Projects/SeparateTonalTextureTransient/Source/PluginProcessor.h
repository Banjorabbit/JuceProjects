/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "SeparateTonalTextureTransient.h"

//==============================================================================
/**
*/
class SeparateTonalTextureTransientAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    SeparateTonalTextureTransientAudioProcessor();
    ~SeparateTonalTextureTransientAudioProcessor();

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

	enum SelectOutput { TONAL, TEXTURE, TRANSIENT, SUM};
	float tonalThreshold;
	void UpdateParameters(const SelectOutput& outputSelector)
	{
		auto p = Separator.GetParameters();
		switch (outputSelector)
		{
		case TONAL: p.OutputSelector = p.TONAL;break;
		case TEXTURE: p.OutputSelector = p.TEXTURE;break;
		case TRANSIENT: p.OutputSelector = p.TRANSIENT;break;
		case SUM: p.OutputSelector = p.SUM;break;
		}
		Separator.SetParameters(p);
		auto pTonal = Separator.TonalSeparator.TonalDetector.GetParameters();
		pTonal.TonalThreshold = tonalThreshold;
		Separator.TonalSeparator.TonalDetector.SetParameters(pTonal);		
	}
private:
	SeparateTonalTextureTransient Separator;
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SeparateTonalTextureTransientAudioProcessor)
};
