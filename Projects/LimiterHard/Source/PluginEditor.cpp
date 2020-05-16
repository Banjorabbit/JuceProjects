/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
LimiterHardAudioProcessorEditor::LimiterHardAudioProcessorEditor (LimiterHardAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);

	PreGain.setSliderStyle(Slider::LinearBarVertical);
	PreGain.setRange(.01f, 10.f, .01f);
	PreGain.setTextValueSuffix(" Pre gain");
	PreGain.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	PreGain.setValue(processor.PreGain);
	addAndMakeVisible(&PreGain);
	PreGain.setComponentID("PreGain");
	PreGain.addListener(this);
	

	PostGain.setSliderStyle(Slider::LinearBarVertical);
	PostGain.setRange(.01f, 10.f, .01f);
	PostGain.setTextValueSuffix(" Post gain");
	PostGain.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	PostGain.setValue(processor.PostGain);
	addAndMakeVisible(&PostGain);
	PostGain.setComponentID("PostGain");
	PostGain.addListener(this);

	LookAhead.setSliderStyle(Slider::LinearBarVertical);
	LookAhead.setRange(.1f, 20.f, .1f);
	LookAhead.setTextValueSuffix(" Look ahead (ms)");
	LookAhead.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	LookAhead.setValue(processor.LookAhead);
	addAndMakeVisible(&LookAhead);
	LookAhead.setComponentID("LookAhead");
	LookAhead.addListener(this);

	Release.setSliderStyle(Slider::LinearBarVertical);
	Release.setRange(.001f, .50f, .001f);
	Release.setTextValueSuffix(" Release (s)");
	Release.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	Release.setValue(processor.Release);
	addAndMakeVisible(&Release);
	Release.setComponentID("Release");
	Release.addListener(this);
}

LimiterHardAudioProcessorEditor::~LimiterHardAudioProcessorEditor()
{
}

//==============================================================================
void LimiterHardAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Limiter!", getLocalBounds(), Justification::centred, 1);
}

void LimiterHardAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	PreGain.setBounds(40, 30, 20, getHeight() - 60);
	PostGain.setBounds(80, 30, 20, getHeight() - 60);
	LookAhead.setBounds(120, 30, 20, getHeight() - 60);
	Release.setBounds(160, 30, 20, getHeight() - 60);
}

void LimiterHardAudioProcessorEditor::sliderValueChanged(Slider * slider)
{
	juce::String name = slider->getComponentID();
	if (name == "PreGain")
	{
		processor.PreGain = static_cast<float>(slider->getValue());
		processor.UpdateParameters();
	}
	else if (name == "PostGain")
	{
		processor.PostGain = static_cast<float>(slider->getValue());
		processor.UpdateParameters();
	}
	else if (name == "LookAhead")
	{
		processor.LookAhead = static_cast<float>(slider->getValue());
		processor.UpdateCoefficients();
	}
	else if (name == "Release")
	{
		processor.Release = static_cast<float>(slider->getValue());
		processor.UpdateParameters();
	}
}
