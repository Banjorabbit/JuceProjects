/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PitchShiftAudioProcessorEditor::PitchShiftAudioProcessorEditor (PitchShiftAudioProcessor& p)
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
}

PitchShiftAudioProcessorEditor::~PitchShiftAudioProcessorEditor()
{
}

//==============================================================================
void PitchShiftAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Pitch Shifter!", getLocalBounds(), Justification::centred, 1);
}

void PitchShiftAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	Pitch.setBounds(40, 30, 20, getHeight() - 60);
}

void PitchShiftAudioProcessorEditor::sliderValueChanged(Slider * slider)
{
	processor.Pitch = static_cast<float>(Pitch.getValue());
	processor.UpdateParameters();
}