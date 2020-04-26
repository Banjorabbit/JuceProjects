/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PitchShiftAdaptiveResolutionPluginAudioProcessorEditor::PitchShiftAdaptiveResolutionPluginAudioProcessorEditor (PitchShiftAdaptiveResolutionPluginAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);

	Pitch.setSliderStyle(Slider::LinearBarVertical);
	Pitch.setRange(.5f, 2.f, .01f);
	Pitch.setTextValueSuffix(" Pitch Factor");
	Pitch.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	Pitch.setValue(processor.Pitch);
	addAndMakeVisible(&Pitch);
	Pitch.addListener(this);

	DetectTransient.setSliderStyle(Slider::LinearBarVertical);
	DetectTransient.setRange(.5f, 1.f, .01f);
	DetectTransient.setTextValueSuffix(" DetectTransient Factor");
	DetectTransient.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	DetectTransient.setValue(processor.DetectTransient);
	addAndMakeVisible(&DetectTransient);
	DetectTransient.addListener(this);
}

PitchShiftAdaptiveResolutionPluginAudioProcessorEditor::~PitchShiftAdaptiveResolutionPluginAudioProcessorEditor()
{
}

//==============================================================================
void PitchShiftAdaptiveResolutionPluginAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Hello World!", getLocalBounds(), Justification::centred, 1);
}

void PitchShiftAdaptiveResolutionPluginAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	Pitch.setBounds(40, 30, 20, getHeight() - 60);
	DetectTransient.setBounds(60, 30, 20, getHeight() - 60);
}

void PitchShiftAdaptiveResolutionPluginAudioProcessorEditor::sliderValueChanged(Slider * slider)
{
	processor.Pitch = static_cast<float>(Pitch.getValue());
	processor.DetectTransient = static_cast<float>(DetectTransient.getValue());
	processor.UpdateParameters();
}