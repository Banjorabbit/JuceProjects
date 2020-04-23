/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
PitchShiftAdaptiveResolutionPluginAudioProcessor::PitchShiftAdaptiveResolutionPluginAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
}

PitchShiftAdaptiveResolutionPluginAudioProcessor::~PitchShiftAdaptiveResolutionPluginAudioProcessor()
{
}

//==============================================================================
const String PitchShiftAdaptiveResolutionPluginAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool PitchShiftAdaptiveResolutionPluginAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool PitchShiftAdaptiveResolutionPluginAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool PitchShiftAdaptiveResolutionPluginAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double PitchShiftAdaptiveResolutionPluginAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int PitchShiftAdaptiveResolutionPluginAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int PitchShiftAdaptiveResolutionPluginAudioProcessor::getCurrentProgram()
{
    return 0;
}

void PitchShiftAdaptiveResolutionPluginAudioProcessor::setCurrentProgram (int index)
{
}

const String PitchShiftAdaptiveResolutionPluginAudioProcessor::getProgramName (int index)
{
    return {};
}

void PitchShiftAdaptiveResolutionPluginAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void PitchShiftAdaptiveResolutionPluginAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
	int nIn = getNumInputChannels();
	PitchShifter0.Initialize(samplesPerBlock, 1, static_cast<float>(sampleRate));
	if (nIn == 2)
	{
		PitchShifter1.Initialize(samplesPerBlock, 1, static_cast<float>(sampleRate));
	}
	setLatencySamples(PitchShifter0.GetLatencyTotalSamples());
}

void PitchShiftAdaptiveResolutionPluginAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool PitchShiftAdaptiveResolutionPluginAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void PitchShiftAdaptiveResolutionPluginAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;

	Eigen::Map<Eigen::ArrayXf> IOMap0(buffer.getWritePointer(0), buffer.getNumSamples());
	PitchShifter0.Process(IOMap0, IOMap0);
	if (getTotalNumOutputChannels() == 2)
	{
		Eigen::Map<Eigen::ArrayXf> IOMap1(buffer.getWritePointer(1), buffer.getNumSamples());
		PitchShifter1.Process(IOMap1, IOMap1);
	}
}

//==============================================================================
bool PitchShiftAdaptiveResolutionPluginAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* PitchShiftAdaptiveResolutionPluginAudioProcessor::createEditor()
{
    return new PitchShiftAdaptiveResolutionPluginAudioProcessorEditor (*this);
}

//==============================================================================
void PitchShiftAdaptiveResolutionPluginAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void PitchShiftAdaptiveResolutionPluginAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new PitchShiftAdaptiveResolutionPluginAudioProcessor();
}
