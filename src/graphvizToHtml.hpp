#pragma once

#include "iosfwd.hpp"
#include "string.hpp"

namespace dinara {

    // This can be used to display a graph in Graphviz format to html.
    // It does the following:
    // - It constructs a Graphviz command to render the graph in svg format.
    // - It invokes the Graphviz command with the given timeout, in seconds.
    // - It writes the svg output to the given html ostream.
    // - Finally, it removes the Graphviz file and the svg file.
    // The options argument should not include the options to specify the
    // input and output file name and the svg format. They should only include
    // rendering options.
    // If successful, this removes the input Graphviz file and the svg file.
    // Otherwise, it throws an exception and does not remove the input graphviz file.
    void graphvizToHtml(
        const string& dotFileName,
        const string& layoutName,
        double timeout,
        const string& options,
        ostream& html);
}

