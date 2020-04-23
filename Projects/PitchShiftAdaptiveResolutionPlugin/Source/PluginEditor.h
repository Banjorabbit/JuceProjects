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
class PitchShiftAdaptiveResolutionPluginAudioProcessorEditor  : public AudioProcessorEditor
{
public:
    PitchShiftAdaptiveResolutionPluginAudioProcessorEditor (PitchShiftAdaptiveResolutionPluginAudioProcessor&);
    ~PitchShiftAdaptiveResolutionPluginAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    PitchShiftAdaptiveResolutionPluginAudioProcessor& processor;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PitchShiftAdaptiveResolutionPluginAudioProcessorEditor)
};
