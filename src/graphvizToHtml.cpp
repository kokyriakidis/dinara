#include "graphvizToHtml.hpp"
#include "runCommandWithTimeout.hpp"
using namespace dinara;

#include <filesystem>
#include "fstream.hpp"
#include "stdexcept.hpp"
#include "string.hpp"

void dinara::graphvizToHtml(
    const string& dotFileName,
    const string& layoutName,
    double timeout,
    const string& options,
    ostream& html)
{
    // Construct the Graphviz command.
    const string svgFileName = dotFileName + ".svg";
    const string command = layoutName + " -T svg " + dotFileName + " -o " + svgFileName + " " + options;

    // Run the command with the given timeout.
    bool timeoutTriggered = false;
    bool signalOccurred = false;
    int returnCode = 0;
    runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
    if(signalOccurred) {
        throw runtime_error("Error during graph layout. Command was " + command);
    }
    if(timeoutTriggered) {
        throw runtime_error("Timeout during graph layout. Command was " + command);
    }
    if(returnCode != 0) {
        throw runtime_error("Error during graph layout. Command was " + command);
    }

    // Success, we can remove the Graphviz file.
    std::filesystem::remove(dotFileName);

    // Write the svg to html.
    html << "<div style='display:inline-block'>";
    ifstream svgFile(svgFileName);
    html << svgFile.rdbuf();
    svgFile.close();
    html << "</div>";

    // Remove the .svg file.
    std::filesystem::remove(svgFileName);
}

