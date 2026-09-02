cmake_minimum_required(VERSION 3.15)
project(DISTORTION_PLUGIN VERSION 0.0.1)

set(CMAKE_CXX_FLAGS_RELWITHASSERT "${CMAKE_CXX_FLAGS_RELEASE} -UNDEBUG")

# 1. Tell CMake where JUCE lives
# Ensure this path points to the root JUCE folder (containing the 'modules' folder)
add_subdirectory(/home/stefan/src/2nd-tree/JUCE _deps/juce)

# 2. Define the Plugin Metadata
juce_add_plugin(DistortionPlugin
    FORMATS VST3
    PRODUCT_NAME "DistortionPlugin"
    COMPANY_NAME "Stefan"
)

target_sources(DistortionPlugin PRIVATE
  DistortionPlugin.cpp
)

target_include_directories(DistortionPlugin PRIVATE /home/stefan/src/pandaresampler/include/)

# 3. Link the Modules and System Flags (THIS is where pkg-config gets triggered)
target_link_libraries(DistortionPlugin PRIVATE
    juce::juce_audio_utils
    juce::juce_audio_processors
    juce::juce_recommended_config_flags
    juce::juce_recommended_lto_flags
    juce::juce_recommended_warning_flags
)

# 4. Fix the VST3 error by passing this to the compiler
target_compile_definitions(DistortionPlugin PRIVATE
    JUCE_VST3_CAN_REPLACE_VST2=0
    JUCE_USE_CURL=0
    JUCE_WEB_BROWSER=0
)
