/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class SeparateTonalTextureTransientAudioProcessorEditor  : public AudioProcessorEditor, private Slider::Listener
{
public:
    SeparateTonalTextureTransientAudioProcessorEditor (SeparateTonalTextureTransientAudioProcessor&);
    ~SeparateTonalTextureTransientAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    SeparateTonalTextureTransientAudioProcessor& processor;

	void sliderValueChanged(Slider *slider) override;
	Slider selectOutput;
	Slider tonalThreshold;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SeparateTonalTextureTransientAudioProcessorEditor)
};
