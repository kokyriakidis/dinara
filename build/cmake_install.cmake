# Install script for directory: /home/kokyriakidis/Downloads/variantClustering/shasta

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ".")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/kokyriakidis/Downloads/variantClustering/shasta/build/wfa2lib/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/kokyriakidis/Downloads/variantClustering/shasta/build/staticLibrary/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/kokyriakidis/Downloads/variantClustering/shasta/build/staticExecutable/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/kokyriakidis/Downloads/variantClustering/shasta/build/dynamicLibrary/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/shasta-install/bin" TYPE PROGRAM FILES
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AlignOrientedReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AlignOrientedReads1.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AlignOrientedReads4.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AnalyzeAlignmentMatrix.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AnalyzeAssemblyGraphBubbles.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AnalyzeReadGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Assemble.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AssembleMarkerGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AssembleMarkerGraphVertices.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/AssembleSegment.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CheckConfigurations.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CheckMarkerGraphIsStrandSymmetric.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CleanupDuplicateMarkers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CleanupRunDirectory.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ClusterMarkerGraphEdgeOrientedReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ColorGfaBySimilarityToSegment.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ColorGfaKeySegments.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ColorGfaWithTwoReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeAlignments.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeAssemblyStatistics.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeMarkerGraphCoverageHistogram.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeMarkerGraphVerticesCoverageData.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeReadGraphConnectedComponents.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ComputeSortedMarkers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Copy.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CopyDirectory.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CountKmers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateAndCleanupMarkerGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateAssembly.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateAssemblyGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateAssemblyGraphVertices.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateCompressedAssemblyGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateConfigurationTable.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateLocalSubgraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerGraphAndTransitiveReduction.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerGraphEdgesStrict.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerGraphSecondaryEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerGraphVertices.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateMarkerKmers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateModules.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreatePrimaryMarkerGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraph2.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraph3.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraph4.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraphMode1.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/CreateReadGraphUsingPseudoPaths.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Detangle.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Detangle2.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FastqGzToFasta.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FastqToFasta.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FastqToFastaAll.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FindAlignmentCandidatesLowHash0.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FindAssemblyGraphBubbles.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FindMarkerGraphReverseComplementEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FindMarkerGraphReverseComplementVertices.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FindMarkers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FlagChimericReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FlagCrossStrandReadGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FlagInconsistentAlignments.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FlagPalindromicReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/FlagPrimaryMarkerGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/GenerateConfig.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/GenerateFeedback.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/GenerateRandomHaplotypes.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/GetConfig.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/GetReadId.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/HistogramReadLength.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/InstallPrerequisites-Ubuntu-Old.sh"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/InstallPrerequisites-Ubuntu.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/InstallPrerequisites-Ubuntu.sh"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Mode2Assembly-A.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Mode2Assembly-B-Prepare.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Mode2Assembly-B.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Mode3Assembly.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/PruneMarkerGraphStrongSubgraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/ReadGraphClustering.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/RemoveReadGraphBridges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/RestoreRun.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/RunAssemblies.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SaveRun.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SetMarkerGraphEdgeFlags.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SetupRunDirectory.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SetupSmallRunDirectory.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SimpleBayesianConsensusCallerCreateBuiltin.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SimplifyMarkerGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SplitMarkerGraphSecondaryEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/StepSequence1.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/SummarizeAssemblies.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/Test.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/TestSimpleBayesianConsensusCaller.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/TransitiveReduction.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/TravisCheckBuildMacOS.sh"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/TravisCheckBuildUbuntu.sh"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteAlignmentCandidates.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteAssemblyGraph.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteBadMarkerGraphVertices.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteFasta.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteGfa.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteGfaBothStrands.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteLocalAlignmentCandidateReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteLocalReadGraphReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteMarkers.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteOrientedRead.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteOrientedReadPath.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteOrientedReadsBySegment.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteParallelMarkerGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WritePseudoPath.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteRead.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteReadGraphEdges.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/WriteReads.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/comparePhaseAssignments.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/dset64Test.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/generateBandageLabelsFromAlignment.py"
    "/home/kokyriakidis/Downloads/variantClustering/shasta/scripts/testGlobalMsa.py"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/shasta-install" TYPE DIRECTORY FILES "/home/kokyriakidis/Downloads/variantClustering/shasta/conf" USE_SOURCE_PERMISSIONS)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/shasta-install" TYPE DIRECTORY FILES "/home/kokyriakidis/Downloads/variantClustering/shasta/docs")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/kokyriakidis/Downloads/variantClustering/shasta/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
