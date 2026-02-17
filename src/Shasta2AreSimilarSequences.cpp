#include "Shasta2AreSimilarSequences.hpp"

#include "Base.hpp"
#include "Shasta2AbpoaWrapper.hpp"
#include "Shasta2IsPeriodic.hpp"
#include "DINARA_ASSERT.hpp"

#include <iterator>

using namespace dinara;
using namespace std;

bool dinara::shasta2AreSimilarSequences(
    const vector<Base>& x,
    const vector<Base>& y,
    const vector<uint64_t>& minRepeatCount,
    ostream& html)
{
    const vector< vector<Base> > sequences = {x, y};
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;

    shasta2Abpoa(sequences, consensus, alignment, alignedConsensus, true);

    DINARA_ASSERT(alignment.size() == 2);
    const vector<AlignedBase>& X = alignment[0];
    const vector<AlignedBase>& Y = alignment[1];
    const uint64_t alignmentLength = X.size();
    DINARA_ASSERT(Y.size() == alignmentLength);

    if(html) {
        html <<
            "<h3>Alignment</h3><div style='font-family:monospace;white-space:nowrap;'>";

        for(size_t position=0; position<alignmentLength; ) {
            if((position % 10) == 0) {
                const string label = to_string(position);
                html << label;
                for(size_t i=0; i<10-label.size(); i++) {
                    html << "&nbsp;";
                }
                position += 10;
            } else {
                html << "&nbsp;";
                ++position;
            }
        }
        html << "<br>";

        for(size_t position=0; position<alignmentLength; position++) {
            if((position % 10) == 0) {
                html << "|";
            } else if((position % 5) == 0) {
                html << "+";
            } else {
                html << ".";
            }
        }

        for(uint64_t i=0; i<2; i++) {
            html << "<br>";
            for(uint64_t position=0; position<alignmentLength; position++) {
                const AlignedBase b = alignment[i][position];
                const AlignedBase bOther = alignment[1 - i][position];
                if(b != bOther) {
                    html << "<span style='background-color:pink'>";
                }
                html << b;
                if(b != bOther) {
                    html << "</span>";
                }
            }
        }
        html << "</div><p>";
    }

    for(uint64_t j=0; j<alignmentLength; j++) {
        const AlignedBase b0 = X[j];
        const AlignedBase b1 = Y[j];
        if((!b0.isGap()) && (!b1.isGap()) && (b0 != b1)) {
            if(html) {
                html << "<br>Mismatches are present, sequences are not similar." << endl;
            }
            return false;
        }
    }

    for(uint64_t j=0; j<2; j++) {
        const vector<AlignedBase>& A = alignment[j];
        const vector<AlignedBase>& B = alignment[1 - j];

        for(uint64_t streakBegin=0; streakBegin<alignmentLength; ) {
            if(!A[streakBegin].isGap()) {
                ++streakBegin;
                continue;
            }
            uint64_t streakEnd = streakBegin + 1;
            for(; streakEnd<alignmentLength; ++streakEnd) {
                if(!A[streakEnd].isGap()) {
                    break;
                }
            }
            const uint64_t streakLength = streakEnd - streakBegin;
            if(html) {
                html << "<p>Found a " << streakLength << " base deletion in the " << (j == 0 ? "first" : "second")
                    << " sequence at alignment positions " << streakBegin << " " << streakEnd << "."
                    << "<br>The deleted sequence is ";
                for(uint64_t position=streakBegin; position!=streakEnd; position++) {
                    html << B[position];
                }
                html << ".";
            }

            for(uint64_t period=1; period<minRepeatCount.size(); period++) {
                if(period > streakLength) {
                    break;
                }
                if(html) {
                    html << "<br>Checking period " << period;
                }
                if(!shasta2IsPeriodic(B.begin()+streakBegin, B.begin()+streakEnd, period)) {
                    if(html) {
                        html << "<br>The deleted sequence does not have period " << period << endl;
                    }
                    continue;
                }

                const uint64_t periodicSequenceBegin = streakBegin;

                uint64_t copyNumberOnRight = 1;
                for(; ; ++copyNumberOnRight) {
                    bool copyIsIntact = true;
                    for(uint64_t i=0; i<period; i++) {
                        const uint64_t shiftedPosition = periodicSequenceBegin + i + period * copyNumberOnRight;
                        if(shiftedPosition >= alignmentLength) {
                            copyIsIntact = false;
                            break;
                        }
                        if(B[shiftedPosition] != B[periodicSequenceBegin + i]) {
                            copyIsIntact = false;
                            break;
                        }
                    }
                    if(!copyIsIntact) {
                        break;
                    }
                }
                --copyNumberOnRight;

                uint64_t copyNumberOnLeft = 1;
                for(; ; ++copyNumberOnLeft) {
                    bool copyIsIntact = true;
                    if(period > streakBegin) {
                        copyIsIntact = false;
                    } else {
                        for(uint64_t i=0; i<period; i++) {
                            const uint64_t shiftedPosition = periodicSequenceBegin + i - period * copyNumberOnLeft;
                            if(B[shiftedPosition] != B[periodicSequenceBegin + i]) {
                                copyIsIntact = false;
                                break;
                            }
                        }
                    }
                    if(!copyIsIntact) {
                        break;
                    }
                }
                --copyNumberOnLeft;

                const uint64_t totalCopyNumber = copyNumberOnRight + copyNumberOnLeft + 1;
                if(html) {
                    html << "<br>Found a total " << totalCopyNumber << " copies of this repeat of period " << period
                        << " (" << totalCopyNumber * period << " bases).";
                }

                if(totalCopyNumber < minRepeatCount[period]) {
                    if(html) {
                        html << "<br>This repeat is short. This difference is significant. These sequences are not similar.";
                    }
                    return false;
                } else {
                    break;
                }
            }

            streakBegin = streakEnd;
        }
    }

    if(html) {
        html << "<br>No significant differences were found. These sequences are similar.";
    }
    return true;
}
