/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
SeparateTonalTextureTransientAudioProcessorEditor::SeparateTonalTextureTransientAudioProcessorEditor (SeparateTonalTextureTransientAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    setSize (400, 300);
	selectOutput.setSliderStyle(Slider::LinearBarVertical);
	selectOutput.setRange(.0f, 3.f, 1.f);
	selectOutput.setTextValueSuffix("Output Selector");
	selectOutput.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	//selectOutput.setValue(processor.Pitch);
	addAndMakeVisible(&selectOutput);
	selectOutput.addListener(this);

	tonalThreshold.setSliderStyle(Slider::LinearBarVertical);
	tonalThreshold.setRange(0.f, 5.f, .1f);
	tonalThreshold.setTextValueSuffix("Tonal Threshold");
	tonalThreshold.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	//selectOutput.setValue(processor.Pitch);
	addAndMakeVisible(&tonalThreshold);
	tonalThreshold.addListener(this);

	timeConstantTransient.setSliderStyle(Slider::LinearBarVertical);
	timeConstantTransient.setRange(0.f, 1.f, .001f);
	timeConstantTransient.setTextValueSuffix("Transient Timeconstant");
	timeConstantTransient.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	//selectOutput.setValue(processor.Pitch);
	addAndMakeVisible(&timeConstantTransient);
	timeConstantTransient.addListener(this);

	timeConstantTonal.setSliderStyle(Slider::LinearBarVertical);
	timeConstantTonal.setRange(0.f, 1.f, .001f);
	timeConstantTonal.setTextValueSuffix("Transient Timeconstant");
	timeConstantTonal.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	//selectOutput.setValue(processor.Pitch);
	addAndMakeVisible(&timeConstantTonal);
	timeConstantTonal.addListener(this);

	predictionDelay.setSliderStyle(Slider::LinearBarVertical);
	predictionDelay.setRange(0.f, 1.f, .001f);
	predictionDelay.setTextValueSuffix("Transient Timeconstant");
	predictionDelay.setTextBoxStyle(Slider::TextBoxBelow, true, 20, getHeight() - 60);
	//selectOutput.setValue(processor.Pitch);
	addAndMakeVisible(&predictionDelay);
	predictionDelay.addListener(this);
}

SeparateTonalTextureTransientAudioProcessorEditor::~SeparateTonalTextureTransientAudioProcessorEditor()
{
}

//==============================================================================
void SeparateTonalTextureTransientAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    g.setColour (Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Hello World!", getLocalBounds(), Justification::centred, 1);
}

void SeparateTonalTextureTransientAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	selectOutput.setBounds(40, 30, 20, getHeight() - 60);
	tonalThreshold.setBounds(60, 30, 20, getHeight() - 60);
	timeConstantTransient.setBounds(80, 30, 20, getHeight() - 60);
	timeConstantTonal.setBounds(100, 30, 20, getHeight() - 60);
	predictionDelay.setBounds(120, 30, 20, getHeight() - 60);
}
void SeparateTonalTextureTransientAudioProcessorEditor::sliderValueChanged(Slider * slider)
{
	processor.tonalThreshold = tonalThreshold.getValue();
	processor.timeConstantTransient = timeConstantTransient.getValue();
	processor.timeConstantTonal = timeConstantTonal.getValue();
	processor.predictionDelay = predictionDelay.getValue();
	SeparateTonalTextureTransientAudioProcessor::SelectOutput select = SeparateTonalTextureTransientAudioProcessor::TRANSIENT;
	if (selectOutput.getValue() < 1.f) { select = SeparateTonalTextureTransientAudioProcessor::TONAL; }
	else if (selectOutput.getValue() < 2.f) { select = SeparateTonalTextureTransientAudioProcessor::TEXTURE; }
	else if (selectOutput.getValue() < 3.f) { select = SeparateTonalTextureTransientAudioProcessor::TRANSIENT; }
	else if (selectOutput.getValue() < 4.f) { select = SeparateTonalTextureTransientAudioProcessor::SUM; }
	processor.UpdateParameters(select);
}